import os
os.chdir("/Users/sarasafciu/Desktop/")

sample_names = []
with open("SPN_raw.vcf", "r") as f:
    while True:
        line = f.readline()
        if not line:
            break
        if line.startswith("#CHROM"):
            sample_names = line.split()[10:]
            break

sample_singleton_pos = {}
for sample in sample_names:
    sample_singleton_pos[sample] = {}

only_bases = {"A", "T", "G", "C"}

with open("SPN_raw.vcf", "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        entries = line.strip().split("\t")
        if len(entries) <10:
            continue
        if entries[6] != "PASS":
            continue
        ref = entries[3]
        alt = entries[4]
        if len(ref) != 1 or len(alt) != 1:
            continue
        if ref not in only_bases or alt not in only_bases:
            continue
        
        sample_genotypes = entries[10:]
        carriers = sum(1 for g in sample_genotypes if g == "1")
        if carriers != 1:
            continue
        contig = entries[0]
        pos = int(entries[1])
        for i, g in enumerate(sample_genotypes):
            if g == "1":
                sample = sample_names[i]
                if contig not in sample_singleton_pos[sample]:
                    sample_singleton_pos[sample][contig] = []
                sample_singleton_pos[sample][contig].append(pos)
total = 0
for sample in sample_singleton_pos.values():
    for contig, positions in sample.items():
        total += len(positions)
print("Total singleton positions stored:", total)

recombinant_pos = set()
for sample in sample_singleton_pos:
    for contig, positions in sample_singleton_pos[sample].items():
        positions.sort()
        for i in range(len(positions)):
            window = []
            for j in range(i + 1, len(positions)):
                if positions[j] - positions[i] <= 2000:
                    window.append(positions[j])
                else:
                    break
            if len(window) >= 3:
                recombinant_pos.add((contig, positions[i]))
                for pos in window:
                    recombinant_pos.add((contig, pos))

print("Total recombinant positions:", len(recombinant_pos))
print("Number of contigs with recombinant positions:", len(set(contig for contig, pos in recombinant_pos)))

with open("SPN_recombinant_snps_test.csv", "w") as f:
    f.write("CHROM,POS\n")
    for contig, pos in recombinant_pos:
          f.write(f"{contig}, {pos}\n")

print("Recombinant SNPs savd to ST_recombinant_snps.csv")
with open("SPN_clean_singletons_test.csv", "w") as f:
    f.write("CHROM,POS\n")
    for sample in sample_singleton_pos:
        for contig, positions in sample_singleton_pos[sample].items():
            for pos in positions:
                if (contig, pos) not in recombinant_pos:
                  f.write(f"{contig}, {pos}\n")
print("Clean singletons saved to ST_clean_singletons.cav")


with open("SPN_raw.vcf", "r") as infile, \
     open("SPN_singletons_only_test.vcf", "w") as singletons_out,\
     open("SPN_recombinant_snps_test.vcf", "w") as recomb_out, \
     open("SPN_clean_singletons_test.vcf", "w") as clean_out:

    for line in infile:
        if line.startswith("#"):
            singletons_out.write(line)
            recomb_out.write(line)
            clean_out.write(line)
            continue
        entries = line.strip().split("\t")
        if len(entries) < 10:
            continue
        if entries[6] != "PASS":
            continue
        ref = entries[3]
        alt = entries[4]
        if len(ref) != 1 or len(alt) != 1:
            continue
        if ref not in only_bases or alt not in only_bases:
            continue

        sample_genotypes = entries[10:]
        carriers = sum(1 for g in sample_genotypes if g == "1")
        if carriers != 1:
            continue
        contig = entries[0]
        pos = int(entries[1])

        singletons_out.write(line)

        if (contig,pos) in recombinant_pos:
            recomb_out.write(line)
        else:
            clean_out.write(line)

print("VCF files saved!")


