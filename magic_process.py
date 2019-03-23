import magic
import scprep
import os
import numpy as np
import pandas as pd
# import matplotlib
# import matplotlib.pyplot as plt

# Matplotlib command for Jupyter notebooks only
# %matplotlib inline

# Tutorial for MAGIC:
# https://nbviewer.jupyter.org/github/KrishnaswamyLab/MAGIC/blob/master/python/tutorial_notebooks/emt_tutorial.ipynb

# txt to csv
os.system("cat GSM2544738_K562_dCas9-KRAB_10_sgRNAs_targeting_TAD2_SE2_enhancers_low_MOI_batch_1.counts.txt | tr -s '[:blank:]' ',' > GSM2544738_K562_dCas9-KRAB_10_sgRNAs_targeting_TAD2_SE2_enhancers_low_MOI_batch_1.counts.csv")

# read data, transpose
emt_data = scprep.io.load_csv('GSM2544738_K562_dCas9-KRAB_10_sgRNAs_targeting_TAD2_SE2_enhancers_low_MOI_batch_1.counts.csv').T

# Drop unnecessary information
emt_data = emt_data.drop(["Chr", "Start", "End", "Strand", "Length"], axis=0)
emt_data = emt_data.drop(["barcode_region;lentiGuide-MS2-puro-barcode-empty;lentiGuide-MS2-puro-barcode-empty;+;568;5488"], axis=1)
emt_data = emt_data.reset_index().drop(['index'], axis=1)

# change type from Object to float
emt_data = emt_data.astype("float64")

# Preprocessing: Filtering
# See https://nbviewer.jupyter.org/github/KrishnaswamyLab/MAGIC/blob/master/python/tutorial_notebooks/emt_tutorial.ipynb
# scprep.plot.plot_library_size(emt_data, cutoff=1500)
emt_data = scprep.filter.filter_library_size(emt_data, cutoff=1500)

# Preprocessing: Normalize
# See https://nbviewer.jupyter.org/github/KrishnaswamyLab/MAGIC/blob/master/python/tutorial_notebooks/emt_tutorial.ipynb
emt_data = scprep.normalize.library_size_normalize(emt_data)
emt_data = scprep.transform.sqrt(emt_data)

# Running Magic
magic_op = magic.MAGIC()
# Choosing the genes (arbitrary number of genes)
gene1 = 'CLIC4;NM_013943;chr1;+;25071759;25170815'
gene2 = 'KIR2DS2;NM_001291696;chr19_gl000209_random;+;131432;145743'
gene3 = 'KIR2DS2;NM_001291701;chr19_gl000209_random;+;131432;145743'
emt_magic = magic_op.fit_transform(emt_data, genes=[gene1, gene2, gene3])


# output
output_file = "magic_output.csv"
emt_magic.to_csv(output_file, encoding='utf-8', index=False)

# # Visualizing
# fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16, 6))

# scprep.plot.scatter(x=emt_data[gene1], y=emt_data[gene2], c=emt_data[gene3],  ax=ax1,
#                     xlabel=gene1, ylabel=gene2, legend_title=gene3, title='Before MAGIC')

# scprep.plot.scatter(x=emt_magic[gene1], y=emt_magic[gene2], c=emt_magic[gene3], ax=ax2,
#                     xlabel=gene1, ylabel=gene2, legend_title=gene3, title='After MAGIC')

# plt.tight_layout()
# plt.show()
