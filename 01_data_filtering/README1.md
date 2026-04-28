The following script:
- Filters for "PASS" SNPs whilst skipping SNPs with multi-allelic sites,
- Isolates singletons (i.e. SNPs which are only found in one genome out of the 1000 genome sample),
- Filters for recombination using a 2000bp window and a 3 SNP threshold, \
- Produces three separate VCF files containing the isolated unique SNPs only,
the isolated unique SNPs flagged for recombination,and the 'clean' recombination-free unique SNPs,

The script I uploaded uses my Salmonella typhimurium vcf file and reference fasta file as an example, however
I also ran this same script on the remaining three bacterial species, producing the same results
