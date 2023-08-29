#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Data::Dumper;
############# GetOptions #################
my ($invcf, $inphy, $outdir, $shdir, $key, $bootnum, $fraction, $model);
my ($goalign, $gotree, $fastme, $plink, $ped2fasta, $vcfutils2, $bootstrap_seq);
my ($config_file, $pattern);
my $run;
GetOptions(
	"h|?"=>\&help,
	"v:s"=>\$invcf,
	"p:s"=>\$inphy,
	"o:s"=>\$outdir,
	"k:s"=>\$key,
	"s:s"=>\$shdir,
	"b:s"=>\$bootnum,
	"f:s"=>\$fraction,
	"m:s"=>\$model,
	"c:s"=>\$config_file,
	"p:s"=>\$pattern,
	"r:s"=>\$run,
) || &help;
&help unless (($invcf || $inphy) && $outdir && $key);

sub help
{
	print"
	Description: vcf to NJ-tree by pre-computed distance matrix and parallelizing bootstraps

	-v  <file>  input vcf file                 [force if no -p]
	-p  <file>  input phylip file              [force if no -v]
	-o  <dir>   output directory               [force]
	-k  <str>   output prefix                  [force]
	-s  <dir>   shell directory                [-o]
	-b  <int>   bootstrap times                [100]
	-f  <float> fraction of bootstrap sampling [0.25]
	-m  <str>   substitution model             [pdist]
                pdist: p-distance
                jc   : Juke-Cantor
                k2p  : Kimura 2 Parameters
                f81  : Felsenstein 81
                f84  : Felsenstein 84
                tn93 : Tamura and Nei 1993
	-c  <file>  config file                    [optional]
	-p  <str>   config split pattern           [:=]
	-r  <T/F>   run or not                     [T]
	-h          Help document
";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
if(defined $invcf){
	$invcf = &abs_dir($invcf);
}
if(defined $inphy){
	$inphy = &abs_dir($inphy);
}
if(!-d $outdir){
	`mkdir -p $outdir`;
}
$outdir = &abs_dir($outdir);
my $tmpdir = "$outdir/tmp";
`mkdir -p $tmpdir`;

$shdir ||= $outdir;
if(!-d $shdir){
	`mkdir -p $shdir`;
}
$shdir = &abs_dir($shdir);

if($config_file){
	$config_file = &abs_dir($config_file);
}

$run ||= "T";
$pattern ||= ":=";
##### read config
my %config;
&read_config($config_file);

$config{snp_tree_bootnum} //= 100;
$config{snp_tree_model} ||= "pdist";
$config{snp_tree_bootfrac} //= 1;
$config{thread_small} ||= 4;
$config{perl} ||= "/Bio/bin/perl";
$config{rscript} ||= "/Bio/User/kxie/software/mambaforge/envs/R/bin/Rscript";
$config{plink} ||= "/Bio/Bin/pipeline/gwas/v2.0/bin/plink2";
$config{goalign} ||= "/Bio/User/kxie/software/mambaforge/envs/goalign/bin/goalign";
$config{gotree} ||= "/Bio/User/kxie/software/mambaforge/envs/goalign/bin/gotree";
$config{fastme} ||= "/Bio/User/kxie/software/mambaforge/envs/fastme/bin/fastme";
$config{maxjob} ||= 20;
$config{queue} ||= "all.q,fast.q,centos7";
$config{submit} ||= "qsub";


&main();

############# end time ###################
$current_T = &date_format(localtime());
print "Programe end: $current_T\n\n";
&Runtime($begin_time);
##########################################
# sub main
sub main()
{
	open(PIPE, ">$shdir/tree.step_by_step.sh");
	print PIPE "set -e\n";
	# filter
	my ($cmd, $phyfile);
	if(defined $invcf && !defined $inphy){
		($cmd, $phyfile) = &filter($invcf);
	}else{
		$phyfile = $inphy;
	}
	print PIPE "$cmd\n";

	# distance
	my ($cmd1, $cmd2) = &distance($phyfile);
	print PIPE "$cmd1\n";
	print PIPE "$cmd2\n";

	# post
	$cmd = &treeformat();
	print PIPE "$cmd\n";
	close(PIPE);

	if($run eq "T"){
		&run_or_die("$shdir/tree.step_by_step.sh");
	}
}
# sub filter
sub filter()
{
	my ($vcf) = @_;
	my $phy = "$outdir/$key.phy.gz";
	open(SH, ">$shdir/tree.q01.filter.sh");
	my $cmd = "$config{perl} $Bin/vcfutils2.pl connect -v $vcf -o $outdir/$key.connect.vcf.gz ";
	$cmd .= "-r $outdir/$key\.connect.rm.vcf.gz -d T -a T -c T -f none";
	$cmd = &srun($config{submit}, $cmd, $config{queue}, "1", $config{thread_small});
	#($slurm, $scmd, $queue, $node, $thread)
	print SH "$cmd\n";

	$cmd = "$config{plink} --vcf $outdir/$key\.connect.vcf.gz --recode ped --snps-only just-acgt ";
	$cmd .= " --double-id --out $outdir/$key\.connect";
	$cmd = &srun($config{submit}, $cmd, $config{queue}, "1", $config{thread_small});
	print SH "$cmd\n";

	$cmd = "$config{perl} $Bin/ped2fasta.pl -p $outdir/$key\.connect.ped -o $phy -f p";
	$cmd = &srun($config{submit}, $cmd, $config{queue}, "1", $config{thread_small});
	print SH "$cmd\n";
	close(SH);
	$cmd = "bash $shdir/tree.q01.filter.sh";
	return($cmd, $phy);
}

sub distance()
{
	my ($phy) = @_;
	my $shq02 = "$shdir/tree.q02.dist.sh";
	my $shq03 = "$shdir/tree.q03.fastme.sh";
	open(SH1, ">$shq02");
	open(SH2, ">$shq03");
	my $cmd = "$config{goalign} compute distance -i $phy -p -t $config{thread_small} ";
	$cmd .= "-m $config{snp_tree_model} -o $outdir/$key.orig.dist";
	$cmd = &srun($config{submit}, $cmd, $config{queue}, "1", $config{thread_small});
	#($slurm, $scmd, $queue, $node, $thread)
	print SH1 "$cmd\n";
	$cmd = "$config{fastme} -i $outdir/$key.orig.dist -o $outdir/$key\.orig.tree";
	print SH2 "$cmd\n";
	for (my $i = 0; $i <= $config{snp_tree_bootnum}; $i++) {
		$i = sprintf("%04d", $i);
		$cmd = "$config{goalign} build distboot -i $phy -p -t 4 -m $config{snp_tree_model} ";
		$cmd .= "-f $config{snp_tree_bootfrac} -o $tmpdir/$key.boot$i.dist";
		$cmd = &srun($config{submit}, $cmd, $config{queue}, "1", $config{thread_small});
		print SH1 "$cmd\n";
		$cmd = "$config{fastme} -i $tmpdir/$key\.boot$i\.dist -o $tmpdir/$key\.boot$i\.tree";
		$cmd = &srun($config{submit}, $cmd, $config{queue}, "1", $config{thread_small});
		print SH2 "$cmd\n";
	}
	close(SH1);
	close(SH2);
	my ($cmd1, $cmd2);
	if($config{submit} eq "qsub"){
		$cmd1 = &qsub($shq02, "1", "5", $config{maxjob}, $config{queue}, "10G", "2");
		$cmd2 = &qsub($shq03, "10", "5", $config{maxjob}, $config{queue}, "5G", "1");
	}else{
		$cmd1 = "$config{parafly} -CPU $config{maxjob} -c $shq02 -failed_cmds $shq02.fail";
		$cmd2 = "$config{parafly} -CPU $config{maxjob} -c $shq03 -failed_cmds $shq03.fail";
	}
	return($cmd1, $cmd2);
}

sub treeformat()
{
	my $rspt = "$shdir/tree.bstformat.r";
	open(RS, ">$rspt");
	my $cmd = "options(scipen=100)\n";
	$cmd .= "library(ape)\n";
	$cmd .= "setwd(\"$outdir\")\n";
	$cmd .= "tr <- read.tree(\"$key.add\.tree\")\n";
	$cmd .= "bst <- as.numeric(tr\$node.label) * 100\n";
	$cmd .= "bst[2:length(bst)] <- round(bst[2:length(bst)], 0)\n";
	$cmd .= "bst[is.na(bst)] <- \"\"\n";
	$cmd .= "tr\$node.label <- bst\n";
	$cmd .= "write.tree(tr, file = \"$outdir/$key.nwk\")";
	print RS "$cmd\n";
	close(RS);

	my $shq04 = "$shdir/tree.q04.format.sh";
	open(SH, ">$shq04");
	$cmd = "cat $tmpdir/$key\.boot*\.tree|grep -v \"\^\$\" > $outdir/$key\.boot.tree";
	print SH "$cmd\n";
	$cmd = "$config{gotree} compute support tbe -b $outdir/$key\.boot.tree ";
	$cmd .= "-i $outdir/$key\.orig.tree -o $outdir/$key\.add.tree";
	print SH "$cmd\n";
	$cmd = "$config{rscript} $rspt";
	print SH "$cmd\n";
	$cmd = "$config{perl} -ne \'s\/[\\d\\.]+e-\\d+\/\${[sprintf \"%.8f\" , \$&]}[0]\/g; print \$_\' ";
	$cmd .= "$outdir/$key\.nwk > $outdir/$key\.nwk1";
	print SH "$cmd\n";
	$cmd = "mv $outdir/$key\.nwk1 $outdir/$key\.nwk";
	print SH "$cmd\n";
	close(SH);
	$cmd = "bash $shq04 > $shq04.log";
	return($cmd);
}

# sub read config
sub read_config()
{
	## reading softwares' path
	open CF, "$config_file" || die $!;
	while(<CF>) {
		chomp;
		next if(/^#|^\s*$/);
		my @x = split /$pattern/;
		$x[0] =~ s/\s+//g;
		if(defined $x[1]) {
			$x[1] =~ s/^\s+//g;
			$x[1] =~ s/\s+$//g;
			$config{$x[0]} = $x[1];
		} else {
			&show_log("WARNNING!! line $. ==> $x[0] <== is undefined!");
		}
	}
	close CF;
}
# sub date format
sub date_format()
{
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $day, $hour, $min, $sec);
}

# sub Runtime
sub Runtime()
{
	my ($begin_time) = @_;
	my $now_time = time();
	my $total_time = $now_time - $begin_time;
	my $sec = 0; my $minu = 0; my $hour = 0;
	if($total_time >= 3600){
		$hour = int($total_time/3600);
		my $left = $total_time % 3600;
		if($left >= 60){
			$minu = int($left/60);
			$sec = $left % 60;
			$total_time = $hour."h\-".$minu."m\-".$sec."s";
		} else {
			$minu = 0;
			$sec = $left;
			$total_time = $hour."h\-".$minu."m\-".$sec."s";
		}
	} else {
		if($total_time >= 60){
			$minu = int($total_time/60);
			$sec = $total_time % 60;
			$total_time = $minu."m\-".$sec."s";
		} else {
			$sec = $total_time;
			$total_time = $sec."s";
		}
	}
	print "Total elapsed time [$total_time]\n\n";
}

# sub absolutely directory
sub abs_dir()
{
	my ($in) = @_;
	my $current_dir = `pwd`;
	chomp($current_dir);
	my $return_dir = "";
	if(-f $in){
		my $in_dir = dirname($in);
		my $in_file = basename($in);
		chdir $in_dir;
		$in_dir = `pwd`;
		chomp($in_dir);
		$return_dir = "$in_dir/$in_file";
	} elsif(-d $in) {
		chdir $in;
		my $in_dir = `pwd`;
		chomp($in_dir);
		$return_dir = $in_dir;
	} else {
		die("ERROR: there is no file or dir called [$in], please check!\n");
	}
	chdir $current_dir;
	return $return_dir;
}

# show log
sub show_log()
{
	my ($text) = @_;
	my $current_time = &date_format(localtime());
	print "$current_time: $text\n";
}

# run or die
sub run_or_die()
{
	my ($cmd) = @_;
	&show_log($cmd);
	my $flag = system($cmd);
	if($flag != 0){
		&show_log("ERROR: command $cmd");
		exit(1);
	}
	&show_log("done\n");
}

# qsub
sub qsub()
{
	my ($sh, $lines, $interval, $maxjob, $queue, $vf, $slots) = @_;
	$maxjob ||= 20;
	$queue ||= "all.q,fast.q,centos7";
	$vf ||= "5G";
	$slots ||= 2;
	$interval ||= 3;
	my $cmd_qsub = "/Bio/bin/qsub-sge.pl --convert no --interval $interval --queue $queue ";
	$cmd_qsub .= "--maxjob $maxjob --lines $lines --slots $slots --resource vf=$vf --reqsub $sh";
	return($cmd_qsub);
}

sub srun()
{
	my ($slurm, $scmd, $queue, $node, $thread) = @_;
	my $rcmd;
	if($slurm eq "srun"){
		$rcmd = "$config{srun} -p $queue -N $node -c $thread $scmd";
	}else{
		$rcmd = $scmd;
	}
	return($rcmd);
}

