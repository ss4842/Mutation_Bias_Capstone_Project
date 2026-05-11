import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

species = ["S. typhimurium", "S. pneumoniae", "S. pyogenes", "V. cholerae"]
colours = ["#f4978e" , "#93e1d8", "#758bfd", "#a8c256"]

ts_counts = [2201,5395, 13316, 17555]
tv_counts = [1535, 2001, 3480, 8232]

fig = plt.figure(figsize=(8, 8))

ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

x = np.array([0, 0.15])

ax1.bar(x,[ts_counts[0], tv_counts[0]], color = colours[0], width=0.1)
ax1.set_xticks(x)
ax1.set_xticklabels(["Transitions", "Transversions"])
ax1.set_title("S. typhimurium", fontsize = 12)
ax1.set_ylabel("Count", fontsize =10)

ax2.bar(x, [ts_counts[1], tv_counts[1]], color = colours[1], width=0.1)
ax2.set_xticks(x)
ax2.set_xticklabels(["Transitions", "Transversions"])
ax2.set_title("S. pneumoniae", fontsize = 12)
ax2.set_ylabel("Count", fontsize = 10)

ax3.bar(x, [ts_counts[2], tv_counts[2]], color = colours[2], width=0.1)
ax3.set_xticks(x)
ax3.set_xticklabels(["Transitions", "Transversions"])
ax3.set_title("S. pyogenes", fontsize = 12)
ax3.set_ylabel("Count", fontsize = 10)

ax4.bar(x, [ts_counts[3], tv_counts[3]], color = colours[3], width=0.1)
ax4.set_xticks(x)
ax4.set_xticklabels(["Transitions", "Transversions"])
ax4.set_title("V. cholerae", fontsize = 12)
ax4.set_ylabel("Count", fontsize = 10)

plt.suptitle("Transition and Transversion Analysis Across Four Bacterial Species", fontsize = 12)
plt.subplots_adjust(hspace=0.35, wspace=0.3, top=0.92)
plt.savefig("Ts_Tv_analysis.png", dpi=300)





