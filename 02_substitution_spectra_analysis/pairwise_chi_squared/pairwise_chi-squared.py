from scipy.stats import chi2_contingency
from itertools import combinations

#6-class counts for each species
ST = [1630, 243, 643, 771, 371, 278]
SPN = [4007, 244, 929, 1388, 511, 317]
SPY = [8728, 456, 1429, 4588, 927, 668]
VC = [12766, 1156, 3946, 4789, 1572, 1558]

species_counts = {
    "ST": ST,
    "SPN": SPN,
    "SPY": SPY,
    "VC": VC
}

print("\nPairwise chi-squared tests:")
print("-" * 40)

for (sp1, counts1), (sp2, counts2) in combinations(species_counts.items(), 2):
    table = [counts1, counts2]
    chi2, p, dof, expected = chi2_contingency(table)
    significant = "SIGNIFICANT" if p < 0.05 else "NOT SIGNIFICANT"
    print(f"{sp1} vs {sp2}: chi2={chi2:.2f}, p={p:.6f} ({significant})")

print("-" * 40)
