import os
os.chdir("/Users/ss4842/Desktop/")

vcf_files = [
    ("ST_raw.vcf", "Original VCF"),
    ("ST_singletons_only.vcf", "Singletons VCF"),
    ("ST_recombinant_snps.vcf", "Recombinant VCF"),
    ("ST_clean_singletons.vcf", "Final VCF")
]

with open("ST_vcf_snp_counts.csv", "w") as out:
    out.write("vcf_file,description,snp_count\n")
    for filename, description in vcf_files:
        snp_count = 0
        only_bases = {"C", "G", "T", "A"}
        with open(filename, "r") as f:
            while True:
                line = f.readline()      
                if not line:
                    break
                if line.startswith("#") or len(line.strip()) == 0:
                    continue
                entries = line.split()
        
                filter_status = entries[6]
                if filter_status != "PASS":
                    continue
                ref = entries[3]
                alt = entries[4]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                if ref not in only_bases or alt not in only_bases:
                    continue
                snp_count += 1
        out.write(f"{filename},{description},{snp_count}\n")

print("Table saved to ST_vcf_snp_counts.csv")

