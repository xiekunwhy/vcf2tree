# vcf2tree
construction NJ-tree from vcf file

# Dependences
perl https://www.perl.org
R with ape package installed
plink2 https://www.cog-genomics.org/plink/2.0
parafly https://github.com/ParaFly/ParaFly
goalign https://github.com/evolbioinfo/goalign
gotree https://github.com/evolbioinfo/gotree
fastme http://www.atgc-montpellier.fr/fastme/binaries.php
All can be intall via conda or mambaforge, mambaforge(https://github.com/conda-forge/miniforge, choose Latest installers with Mamba in the base environment:) is recommended for China mainland users.

# example
perl vcf2tree.pl -v Xian.vcf.gz -o vcf2tree -k snp -r F -c software.config.txt
than run q01 to q04 step by step, please parallize q02 to speed up the pipeline.
Change parameters and software path in software.config.txt according to your environment.
