#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
use Data::Dumper;

############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
&main;
############# end time ###################
$current_T = &date_format(localtime());
print "Programe end: $current_T\n\n";
&Runtime($begin_time);
##########################################
exit;

#
sub main {
	&usage if (@ARGV < 1);
	my $command = shift(@ARGV);
	my %func = (subsam=>\&subsam, subsites=>\&subsites, listsam=>\&listsam, rename=>\&rename, 
		connect=>\&connect, popfilter=>\&popfilter);
	die("Unknown command \"$command\".\n") if (!defined($func{$command}));
	&{$func{$command}};
}

sub usage {
	die(qq/
Usage:   vcfutils2.pl <command> [<arguments>]\n
Command: subsam       get a subset of samples
         subsites     get a subset of variant sites
         listsam      list and print sample names
         rename       rename samples
         connect      connect all sites into one fake chromosome
         popfilter    filter vcf in population levels
\n/);
}

sub subsam {
	my ($infile, $outfile, $samplelist);
	GetOptions("v:s"=>\$infile, "o:s"=>\$outfile, "l:s"=>\$samplelist);
	die("Usage: vcfutils2.pl subsam [options]

	-v   <file>   input vcf file      [force]
	-o   <file>   output vcf file     [force]
	-l   <file>   sample list         [force]
\n") unless ($infile && $outfile && $samplelist);
	my %keep;
	open(IN, $samplelist);
	while (<IN>) {
		chomp;
		my @l = split(/\s+/, $_);
		$keep{$l[0]} = 1;
	}
	close(IN);
	my ($fi, $fo);
	open($fi, ($infile =~ /\.gz$/) ? "gzip -dc $infile|" : $infile) || die;
	open($fo, ($outfile =~ /\.gz$/) ? "|gzip > $outfile" : ">$outfile") || die;
	my @col;
	while (<$fi>) {
		chomp;
		if(/^##/){
			print $fo "$_\n";
			next;
		}
		my @t = split(/\t/, $_);
		if(/^#CHROM/){
			for (my $i = 0; $i < @t; $i++) {
				if($i <= 8){
					push(@col, $i);
				}else{
					if($keep{$t[$i]}){
						push(@col, $i);
					}
				}
			}
		}
		print $fo join("\t", map {$t[$_]} @col), "\n";
	}
	close($fi);
	close($fo);
}

sub subsites{
	my ($infile, $outfile, $sitelist, $type);
	GetOptions("v:s"=>\$infile, "o:s"=>\$outfile, "l:s"=>\$sitelist, "t:s"=>\$type);
	die("Usage: vcfutils2.pl subsites [options]

	-v   <file>   input vcf file      [force]
	-o   <file>   output vcf file     [force]
	-l   <file>   site list(pos/id)   [force]
	-t   <pos/id> pos or id           [pos]
\n") unless ($infile && $outfile && $sitelist);
	$type ||= "pos";
	my %keep;
	open(IN, $sitelist);
	while (<IN>) {
		chomp;
		my @l = split(/\s+/, $_);
		my $id = "$l[0]";
		if($type eq "pos"){
			$id .= "__$l[1]";
		}
		$keep{$id} = 1;
	}
	close(IN);

	my ($fi, $fo);
	open($fi, ($infile =~ /\.gz$/)? "gzip -dc $infile |" : $infile) || die;
	open($fo, ($outfile =~ /\.gz$/)? "|gzip > $outfile" : ">$outfile") || die;
	while (<$fi>) {
		chomp;
		if(/^#/){
			print $fo "$_\n";
			next;
		}
		my @l = split(/\t/, $_, 4);
		my $id;
		if($type eq "pos"){
			$id = "$l[0]\__$l[1]";
		}else{
			$id = $l[2];
		}
		if(exists $keep{$id}){
			print $fo "$_\n";
		}
	}
	close($fi);
	close($fo);
}

sub listsam {
	die(qq/Usage: vcfutils2.pl listsam <in.vcf> <out.txt>\n/) if (@ARGV != 2 && -t STDIN);
	my $fn = shift(@ARGV);
	my $fo = shift(@ARGV);
	my $fi;
	open($fi, ($fn =~ /\.gz$/)? "gzip -dc $fn |" : $fn) || die;
	open(OUT, ">$fo");
	while (<$fi>) {
		if (/^#/ && !/^##/) {
			my @t = split;
			for (my $i = 9; $i < @t; $i++) {
				my $j = $i + 1;
				print "$i\t$t[$i]\n";
				print OUT "$i\t$t[$i]\n";
			}
			last;
		}
	}
	close($fi);
	close(OUT);
}

sub rename{
	my ($infile, $outfile, $idmap);
	GetOptions("v:s"=>\$infile, "o:s"=>\$outfile, "d:s"=>\$idmap);
	die("Usage: vcfutils2.pl rename [options]

	-v  <file>  input vcf file     [force]
	-o  <file>  output vcf file    [force]
	-m  <file>  idmap(old new)     [force]
\n") unless ($infile && $outfile && $idmap);
	my %idhash;
	open(IN, $idmap);
	while (<IN>) {
		chomp;
		my ($old, $new) = split(/\s+/, $_);
		$idhash{$old} = $new;
	}
	close(IN);
	my ($fi, $fo);
	open($fi, ($infile =~ /\.gz$/) ? "gzip -dc $infile|" : $infile);
	open($fo, ($outfile =~ /\.gz$/) ? "|gzip > $outfile" : ">$outfile");
	while (<IN>) {
		chomp;
		if(/^#CHROM/){
			my @new;
			my @t = split(/\t/, $_);
			push(@new, @t[0..8]);
			for (my $i = 9; $i < @t; $i++) {
				if(!exists $idhash{$t[$i]}){
					die("ERROR: sample $t[$i] not found in idmap file, please check.\n");
				}else{
					push(@new, $idhash{$t[$i]});
				}
			}
			my $newline = join("\t", @new);
			print $fo "$newline\n";
		}
		print $fo "$_\n";
	}
	close($fi);
	close($fo);
}

sub connect {
	my ($infile, $outfile, $rmfile, $dedup, $addid, $connection, $pattern, $flag);
	GetOptions("v:s"=>\$infile, "o:s"=>\$outfile, "r:s"=>\$rmfile,"d:s"=>\$dedup,
		"a:s"=>\$addid, "c:s"=>\$connection, "p:s"=>\$pattern, "f:s"=>\$flag);
	die("Usage: vcfutils2.pl connect [options]

	-v  <file>  input vcf file                 [force]
	-o  <file>  output vcf file                [force]
	-r  <file>  remove sites outfile           [force]
	-d  <T/F>   remove duplication positions   [T]
	-a  <T/F>   add id chr__pos                [T]
	-c  <T/F>   connect all in one chromosome  [T]
	-p  <str>   connect pattern                [__]
	-f  <str>   filter flag to keep            [PASS]
if remove duplication positions, only first variance will be keept
\n") unless ($infile && $outfile && $rmfile);

	$dedup ||= "T";
	$addid ||= "T";
	$connection ||= "T";
	$pattern ||= "__";
	$flag ||= "PASS";
	my ($fi, $fo, $rm);
	open($fi, ($infile =~ /\.gz$/) ? "gzip -dc $infile|" : $infile) || die;
	open($fo, ($outfile =~ /\.gz$/) ? "|gzip > $outfile" : ">$outfile") || die;
	open($rm, ($rmfile =~ /\.gz$/) ? "|gzip > $rmfile" : ">$rmfile") || die;
	my %ids;
	my $n = 1;
	while (<$fi>) {
		chomp;
		if(/^#/){
			print $fo "$_\n";
			print $rm "$_\n";
		}else{
			my @t = split(/\t/, $_, 10);

			if($t[2] eq "\." && $addid eq "T"){
				$t[2] = $t[0] . $pattern . $t[1];
			}elsif($t[2] ne '\.' and $addid eq 'T') {
				$t[2] = $t[0] . $pattern . $t[1] . $pattern . $t[2];
			}
			if($dedup eq "T" && exists $ids{$t[2]}){
				print $rm join("\t", @t), "\n";
				next;
			}
			if($flag ne "none" && $flag ne $t[6]){
				print $rm join("\t", @t), "\n";
				next;
			}
			$ids{$t[2]} = 1;
			if($connection eq "T"){
				$t[0] = 1;
				$t[1] = $n;
			}
			print $fo join("\t", @t), "\n";
			$n++;
		}
	}
	close($fi);
	close($fo);
}

#### popfilter start
my ($infile, $outfile);
my ($max_miss, $min_maf, $min_mac, $max_het, $min_g1_g2_rate, $pv_hwe);
my ($min_qual, $filter_flag, $min_mean_depth, $max_per_batch);
my ($biallelic, $rmIndel, $rmSnp);
my $tot_indiv;
my $threads;
my $sort;
sub popfilter {
	use threads;
	use threads::shared;
	GetOptions("v:s"=>\$infile, "o:s"=>\$outfile, "f:s"=>\$filter_flag, "p:s"=>\$pv_hwe,
		"t:s"=>\$threads, "s:s"=>\$sort, "ba"=>\$biallelic, "ri"=>\$rmIndel, "rs"=>\$rmSnp, "mm:s"=>\$max_miss,
		"ma:s"=>\$min_maf, "mc:s"=>\$min_mac, "mh:s"=>\$max_het, "hh:s"=>\$min_g1_g2_rate,
		"mq:s"=>\$min_qual, "md:s"=>\$min_mean_depth, "pb:s"=>\$max_per_batch
	);
	die("Usage:   vcfutils2.pl popfilter [options]

	-v      <file>      input vcf file             [force]
	-o      <file>      output vcf file            [force]
	-f      <str>       filter flag(PASS)          [undefined]
	-p      <float>     max hwe p value            [undefined]
	-t      <int>       threads number             [4]
	-s      <any>       sort after filtering       [F]
	-ba     <any>       bi-allelic only            [undefined]
	-ri     <any>       remove indels              [undefined]
	-rs     <any>       remove snp                 [undefined]
	-mm     <float>     max missing rate           [0.2]
	-ma     <float>     min maf                    [0.01]
	-mc     <int>       min mac                    [0]
	-mh     <float>     max het rate               [0.8]
	-hh     <float>     min hom1/hom2|hom2/hom1    [undefined]
	-mq     <float>     min qual                   [undefined]
	-md     <float>     min mean depth             [undefined]
	-pb     <int>       max line per batch         [9999]
\n") unless ($infile && $outfile);
	$infile = &abs_dir($infile);
	my $outdir = dirname($outfile);
	$outdir = &abs_dir($outdir);
	my $outbase = basename($outfile);
	$outfile = "$outdir/$outbase";
	$max_miss //= 0.2;
	$min_maf //= 0.01;
	$min_mac //= 0;
	$max_het //= 0.8;
	$max_per_batch ||= 9999;
	$threads ||= 4;
	$sort ||= "F";
	my %out :shared;
	my ($cnt, $snp_total) :shared;
	$cnt = 0;
	$snp_total = 0;
	my ($fi, $fo);
	open($fi, ($infile =~ /\.gz$/) ? "gzip -dc $infile|" : $infile) || die;
	open($fo, ($outfile =~ /\.gz$/) ? "|gzip > $outfile" : ">$outfile") || die;
	while (<$fi>) {
		chomp;
		if(/^#/){
			print $fo "$_\n";
			if(/^#CHROM/){
				my @t = split(/\t/, $_);
				$tot_indiv = scalar(@t) - 9;
			}
			next;
		}
		my @batch = ($_);
		my @self  = ();
		for (1..$threads) {
			for (1..$max_per_batch) {
				last if eof;
				my $vcf = <$fi>;
				push @batch, $vcf; 
			}
			my $t = threads->create('filter', \@batch, \%out, \$cnt, \$snp_total);
			push @self, $t;
			undef @batch;
			last if eof;
		}
		map { $_->join();} @self;
		print "Writing file... $cnt of $snp_total\r" if $cnt;
		my @ids = &nsort(keys %out);
		foreach my $k (@ids) {
			print $fo "$out{$k}\n";
		}
		undef %out;
	}
	close($fi);
	close($fo);
	if ($sort eq "T") {
		print "Sorting file...";
		my $cmd_sort;
		if (substr($outfile, -2) eq 'gz') {
			$cmd_sort = "(zgrep '^#' $outfile; zgrep -v '^#' $outfile | sort -V -T $outdir -k1,1 -k2,2n ) | gzip -c >$outfile.1 && mv $outfile.1 $outfile";
		} else {
			$cmd_sort = "(grep '^#' $outfile; grep -v '^#' $outfile | sort -V -T $outdir -k1,1 -k2,2n) >$outfile.1 && mv $outfile.1 $outfile";
		}
		&run_or_die($cmd_sort);
		print "done.\n";
	}
	print "=================================================================================";
	print "\nFinal retained\t: $cnt SNPs of $tot_indiv individuals\nTotal $snp_total SNPs\n";
}

# sub filter
sub filter {
	my ($vcf, $out, $cnt, $snp_total) = @_;
	foreach my $vcfline (@$vcf) {
		#################################################################
		#VCF format:
		#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Indiv_genotypes ...
		#Genotypes:GT:PL:DP:SP:GQ ...
		
		#Filter order: biallelic -> Q -> depth -> GQ -> coverage -> local Ho or Fis or MAF -> hwe
		#-> missing rate -> global Ho or Fis or MAF.
		
		#################################################################
		chomp $vcfline;
		if (1) {
			lock($$snp_total);
			$$snp_total++; # Number of total SNPs.
		}
		my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/, $vcfline);
		next if ($tot_indiv != @genos);
		next if ($ref =~ /\*/ || $alt =~ /\*/);
		next if ($biallelic && $alt =~ /\,|\./);
		next if ($rmIndel && ($ref =~ /\w{2,}|\-/ || $alt =~ /\w{2,}|\-/));
		next if ($rmSnp && ($ref !~ /\w{2,}|\-/ && $alt !~ /\w{2,}|\-/));
		next if ($min_qual && $qual < $min_qual);
		next if ($filter_flag && $filter ne $filter_flag);
		my ($gti, $adi, $dpi) = &field_format($format);
		my @gts_total        = ();
		my $tot_miss         = 0;
		my $cnt_Ho           = 0;
		my $cnt_hwe          = 0;
		my $tDP              = 0;
		my $recode = 0;
		for (my $i = 0; $i < @genos; $i++) {
			my $geno0 = $genos[$i];
			my @geno = split(/:/, $geno0);
			#### gt
			my $GT = $geno[$gti];
			push (@gts_total, $GT);
			if($GT eq "./." || $GT eq ".|." || $GT eq "\."){
				$tot_miss++;
				next;
			}
			#### dp
			if($min_mean_depth){
				my $DP;
				if(defined $dpi){
					$DP = $geno[$dpi];
				}elsif(defined $adi){
					my @ads = split(/,/, $geno[$adi]);
					$DP = $ads[0] + $ads[1];
				}else{
					$DP = 50;
				}
				$tDP += $DP;
			}
		}
		###### Global missing rate #######
		my $g_N = $tot_indiv - $tot_miss; # global num of individuals.
		my $miss_rate = $tot_miss/$tot_indiv;
		next if ($miss_rate > $max_miss);
		###### total depth and average depth ######
		if (defined $min_mean_depth) {
			my $aDP = $tDP / $g_N; # average depth.
			next if ($aDP < $min_mean_depth);
		}
		###### Global ######
		my ($Ho, $maf, $mac, $pvalue_hwe, $hom1_hom2) = &calc_stat(\@gts_total, $ref, $alt);
		next if ($maf < $min_maf);
		next if ($mac < $min_mac);
		next if ($Ho > $max_het);
		next if (defined $min_g1_g2_rate && $hom1_hom2 < $min_g1_g2_rate);
		next if (defined $pv_hwe && $pvalue_hwe < $pv_hwe);
		###### Print each snp sites ######
		#lock(@out);
		lock(%{$out});
		lock($$cnt);
		#push @out, join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @filtered) . "\n";
		my $id_key = "$chrom\__$pos";
		if($id eq "\."){
			$id = $id_key;
		}
		${$out}{$id_key} = join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos);
		$$cnt++; # SNP number.
	}
}

sub field_format()
{
	my ($fmat) = @_;
	my ($gti, $adi, $dpi);
	my @fm = split(/:/, $fmat);
	for (my $i = 0; $i < @fm; $i++) {
		if($fm[$i] eq "GT") {
			$gti = $i;
			next;
		}
		if($fm[$i] eq "AD"){
			$adi = $i;
		}
		if($fm[$i] eq "DP"){
			$dpi = $i;
		}
	}
	if(!defined $gti){
		&show_log("ERROR: GT is not included in FORMAT, please check!");
		die;
	}
	return ($gti, $adi, $dpi);
}

sub calc_stat {
	my ($gts, $ref, $alt) = @_;
	my @ref_alt;
	push(@ref_alt, $ref);
	if($alt =~ /,/){
		my @alt_split = split(/,/, $alt);
		push(@ref_alt, @alt_split);
	}else{
		push(@ref_alt, $alt);
	}
	my %allele;
	my %genotype;
	for (my $i = 0; $i < @ref_alt; $i++) {
		$allele{$i} = 0;
		for (my $j = $i; $j < @ref_alt; $j++) {
			$genotype{$i . $j} = 0;
		}
	}
	my $n_indv = 0;
	foreach my $gt (@$gts) {
		next if($gt =~ /\./);
		$n_indv++;
		my @tmp = split(/\/|\|/, $gt);
		@tmp = sort(@tmp);
		$allele{$tmp[0]}++;
		$allele{$tmp[1]}++;
		$genotype{$tmp[0] . $tmp[1]}++;
	}
	my @ref_alt_sorted = (sort{$allele{$b} <=> $allele{$a}} keys %allele);
	my ($mac, $maf);
	my ($obs_hets, $obs_hom1, $obs_hom2);
	my $major = $ref_alt_sorted[0];
	$obs_hom1 = $genotype{$major . $major};
	if(@ref_alt_sorted < 2){
		$mac = 0;
		$maf = 0;
		$obs_hets = 0;
		$obs_hom2 = 0;
	}else{
		$mac = $allele{$ref_alt_sorted[1]};
		$maf = $mac/($n_indv * 2);
		my $minor = $ref_alt_sorted[1];
		if($major > $minor){
			$obs_hets = $genotype{$minor . $major};
		}else{
			$obs_hets = $genotype{$major . $minor};
		}
		$obs_hom2 = $genotype{$minor . $minor};
	}
	if(!defined $obs_hom2){
		print Dumper %genotype;
		die;
	}
	my $pvalue_hwe = 1;
	if(defined $pv_hwe){
		$pvalue_hwe = &hwe($obs_hets, $obs_hom1, $obs_hom2);
	}
	my $hom_rate;
	if($obs_hom1 == $obs_hom2 && $obs_hom1 == 0){
		$hom_rate = 0;
	}elsif($obs_hom1 >= $obs_hom2){
		$hom_rate = $obs_hom2/$obs_hom1;
	}elsif($obs_hom1 < $obs_hom2){
		$hom_rate = $obs_hom1/$obs_hom2;
	}
	my $Ho = $obs_hets/($obs_hets + $obs_hom1 + $obs_hom2);
	return($Ho, $maf, $mac, $pvalue_hwe, $hom_rate);
}

sub hwe {
	my $obs_hets = shift;
	my $obs_hom1 = shift;
	my $obs_hom2 = shift;
	if($obs_hom1 < 0 || $obs_hom2 < 0 || $obs_hets <0) {
		return(-1);
	}
	# rare homozygotes
	my $obs_homr;
	# common homozygotes
	my $obs_homc;
	if($obs_hom1 < $obs_hom2) {
		$obs_homr = $obs_hom1;
		$obs_homc = $obs_hom2;
	} else {
		$obs_homr = $obs_hom2;
		$obs_homc = $obs_hom1;
	}
	# number of rare allele copies
	my $rare_copies = 2 * $obs_homr + $obs_hets;
	# total number of genotypes
	my $genotypes = $obs_homr + $obs_homc + $obs_hets;
	if($genotypes <= 0) {
		return(-1);
	}
	# Initialize probability array
	my @het_probs;
	for(my $i=0; $i<=$rare_copies; $i++) {
		$het_probs[$i] = 0.0;
	}
	# start at midpoint
	my $mid = int($rare_copies * (2 * $genotypes - $rare_copies) / (2 * $genotypes));
	# check to ensure that midpoint and rare alleles have same parity
	if(($rare_copies & 1) ^ ($mid & 1)) {
		$mid++;
	}
	my $curr_hets = $mid;
	my $curr_homr = ($rare_copies - $mid) / 2;
	my $curr_homc = $genotypes - $curr_hets - $curr_homr;
	$het_probs[$mid] = 1.0;
	my $sum = $het_probs[$mid];
	for($curr_hets = $mid; $curr_hets > 1; $curr_hets -= 2) {
		$het_probs[$curr_hets - 2] = $het_probs[$curr_hets] * $curr_hets * ($curr_hets - 1.0) / (4.0 * ($curr_homr + 1.0) * ($curr_homc + 1.0));
		$sum += $het_probs[$curr_hets - 2];
		# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
		$curr_homr++;
		$curr_homc++;
	}
	$curr_hets = $mid;
	$curr_homr = ($rare_copies - $mid) / 2;
	$curr_homc = $genotypes - $curr_hets - $curr_homr;
	for($curr_hets = $mid; $curr_hets <= $rare_copies - 2; $curr_hets += 2) {
		$het_probs[$curr_hets + 2] = $het_probs[$curr_hets] * 4.0 * $curr_homr * $curr_homc / (($curr_hets + 2.0) * ($curr_hets + 1.0));
		$sum += $het_probs[$curr_hets + 2];

		# add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
		$curr_homr--;
		$curr_homc--;
	}
	for(my $i=0; $i<=$rare_copies; $i++) {
		$het_probs[$i] /= $sum;
	}
	# Initialise P-value 
	my $p_hwe = 0.0;
	# P-value calculation for p_hwe
	for(my $i = 0; $i <= $rare_copies; $i++) {
		if($het_probs[$i] > $het_probs[$obs_hets]) {
			next;
		}
		$p_hwe += $het_probs[$i];
	}
	if($p_hwe > 1) {
		$p_hwe = 1.0;
	}
	return($p_hwe);
}

#### popfilter end

sub nsort {
	return @_ if @_ < 2;   # Just to be CLEVER.
	my($x, $i);  # scratch vars
	map
		$_->[0],
	sort {
		# Uses $i as the index variable, $x as the result.
		$x = 0;
		$i = 1;
		while($i < @$a and $i < @$b) {
			last if ($x = ($a->[$i] cmp $b->[$i])); # lexicographic
			++$i;
			last if ($x = ($a->[$i] <=> $b->[$i])); # numeric
			++$i;
		}
		$x || (@$a <=> @$b) || ($a->[0] cmp $b->[0]);
	}
	map {
		my @bit = ($x = defined($_) ? $_ : '');
		if($x =~ m/^[+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee](?:[+-]?\d+))?\z/s) {
			# It's entirely purely numeric, so treat it specially:
			push @bit, '', $x;
		} else {
			# Consume the string.
			while(length $x) {
				push @bit, ($x =~ s/^(\D+)//s) ? lc($1) : '';
				push @bit, ($x =~ s/^(\d+)//s) ?    $1  :  0;
			}
		}
		\@bit;
	}
	@_;
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
