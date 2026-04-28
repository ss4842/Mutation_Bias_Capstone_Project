Once again, all scripts uploaded use Salmonella typhimurium as an example, although I did use the same exact script on all four of my allocated bacterial species.

In the "substitution_analysis.py" file:
- I used my the vcf file in which I stored recombination-free singletons (XX_clean_singletons.vcf)
- from this, I counted each of the 12 possible substitution types per genome, producing the "XX_substitution_counts.csv", which would hold all genome fasta codes (e.g. SAMN30941302.fa) in the first column, and each of the 12 possible substitution counts in the 12 remaining columns, each row holding the total substitution counts for each particular genome.
- I then produced a substitution count matrix ("XX_substitution_count_matrix.csv)
- To then normalise my substitution counts, accounting for their GC and AT genome content,
  - I calculated the GC and AT content of the reference genome by counting the number     of each valid base, then adding together the number of G+C and A+T bases
  - I accounted for the double-stranded nature of DNA, adding together the
    counts of substitutions which would appear as pairs (e.g. C>T/G>A) on DNA (producing 6-class           substitution spectra
  - I divided each C:G substitution by the total GC content and each T:A substitution     by the total AT content in the reference genome for that species,
