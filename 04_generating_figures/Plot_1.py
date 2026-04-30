import os
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/Users/sarasafciu/Desktop/csvs/")

st = pd.read_csv("ST_sub_counts_plot.csv")
spn = pd.read_csv("SPN_sub_counts_plot.csv")
spy = pd.read_csv("SPY_sub_counts_plot.csv")
vc = pd.read_csv("VC_sub_counts_plot.csv")
                  
fig, axes = plt.subplots(2, 2, figsize =(14, 10))

axes[0, 0].bar(st["substitution"], st["count"], color= "#f4978e")
axes[0, 0].set_title("S. typhimurium", fontsize = 14)

axes[0, 1].bar(spn["substitution"], spn["count"], color = "#93e1d8")
axes[0, 1].set_title("S. pneumoniae", fontsize = 14)

axes[1, 1].bar(spy["substitution"], spy["count"], color = "#758bfd")
axes[1, 1].set_title("S. pyogenes", fontsize = 14)

axes[1, 0].bar(vc["substitution"], vc["count"], color = "#a8c256")
axes[1, 0].set_title("V. cholerae", fontsize = 14)


for ax in axes.flat:
    ax.set_xlabel("Substitution Type", fontsize = 12)
    ax.set_ylabel("Mutation Frequency", fontsize = 12)
    ax.tick_params(axis="x", rotation = 45, labelsize = 10)

plt.suptitle("12-Class Substitution Spectrum Across Four Bacterial Species", fontsize = 16)
plt.tight_layout()
plt.savefig("substitution_spectrum_4_species_2.png", dpi=300)
plt.show()
