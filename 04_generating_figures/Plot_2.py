import os
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

os.chdir("/Users/sarasafciu/Desktop/")

renormal = pd.read_csv("6_class_sub_spectra_renormalised_plot.csv")
permilbas = pd.read_csv("6_class_sub_spectra_per_mil_bases_plot.csv")
print(renormal.head())
print(permilbas.head())

species = ["ST", "SPN", "SPY", "VC"]
colours = ["#f4978e", "#93e1d8", "#758bfd", "#a8c256"]
labels = ["S. typhimurium", "S. pneumoniae", "S. pyogenes", "V. cholerae"]

x= np.arange(len(renormal["substitution"]))
width = 0.2

fig, axes = plt.subplots(1, 2, figsize = (16,6))

for i, (sp, col, lab) in enumerate(zip(species, colours, labels)):
    axes[0].bar(x+i*width, renormal[sp], width, label = lab, color=col)
    axes[0].set_xticks(x+width*1.5)
    axes[0].set_xticklabels(renormal["substitution"], rotation = 45, fontsize=10)
    axes[0].set_ylabel("Normalised Mutation Frequency", fontsize = 12)
    axes[0].set_title("Opportunity-normalised Substitution Spectrum", fontsize = 13)
    axes[0].legend(fontsize=9)

for i, (sp, col, lab) in enumerate(zip(species, colours, labels)):
    axes[1].bar(x+i*width, permilbas[sp], width, label=lab, color=col)
    axes[1].set_xticks(x+width*1.5)
    axes[1].set_xticklabels(permilbas["substitution"], rotation = 45, fontsize = 10)
    axes[1].set_ylabel("Normalised Mutations per Million Bases", fontsize = 12)
    axes[1].set_title("Opportunity-normalised Substitution Spectrum (per million bases)", fontsize=13)
    axes[1].legend(fontsize = 9)

    plt.suptitle("6-Class Substitution Spectrum Across Four Bacterial Species", fontsize = 15)
    plt.subplots_adjust(top=0.88)
    plt.tight_layout()
    plt.savefig("6_class_substitution_comparison.png", dpi=300)
    plt.show()
    
