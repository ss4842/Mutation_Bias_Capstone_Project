import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
os.chdir("/Users/ss4842/Desktop")

def get_sample_counts(vcf_file):
    sample_names = []
    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#CHROM"):
                sample_names = line.split()[9:]
                break

    counts = [0] * len(sample_names)
    only_bases = {"A", "T", "C", "G"}

    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#") or len(line.strip()) == 0:
                continue
            entries = line.split()
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
            for i, sample in enumerate(entries[9:]):
                if sample == "1":
                    counts[i] += 1

    return counts

st_counts = get_sample_counts("ST_clean_singletons.vcf")
spn_counts = get_sample_counts("SPN_clean_singletons.vcf")
spy_counts = get_sample_counts("SPY_clean_singletons.vcf")
vc_counts = get_sample_counts("VC_clean_singletons.vcf")

print("ST samples:", len(st_counts), "Max:", max(st_counts) if st_counts else "EMPTY")
print("SPN samples:", len(spn_counts), "Max:", max(spn_counts) if spn_counts else "EMPTY")
print("SPY samples:", len(spy_counts), "Max:", max(spy_counts) if spy_counts else "EMPTY")
print("VC samples:", len(vc_counts), "Max:", max(vc_counts) if vc_counts else "EMPTY")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

colours = ["#f4978e", "#93e1d8", "#758bfd", "#a8c256"]
species = ["S. typhimurium", "S. pneumoniae", "S. pyogenes", "V. cholerae"]

panels = [
    (axes[0, 0], st_counts, colours[0], species[0]),
    (axes[0, 1], spn_counts, colours[1], species[1]),
    (axes[1, 0], spy_counts, colours[2], species[2]),
    (axes[1, 1], vc_counts, colours[3], species[3])
]

for ax, counts, colour, name in panels:
    ax.scatter(range(len(counts)), sorted(counts), color=colour, s=10)
    ax.set_title(name, fontsize=12)
    ax.set_xlabel("Sample (ranked by mutation count)", fontsize=10)
    ax.set_ylabel("Number of SNPs", fontsize=10)

plt.suptitle("Per-sample SNP counts across four bacterial species", fontsize=14)
plt.subplots_adjust(hspace=0.4, wspace=0.3, top=0.92)
plt.savefig("per_sample_snp_counts.png", dpi=300)
print("Figure saved!")
