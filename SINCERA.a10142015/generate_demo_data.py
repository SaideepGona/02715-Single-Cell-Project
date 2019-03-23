'''
Author: Saideep Gona

In order to use the mosaic-seq expression data with SINCERA's pipeline demo
script, several datasets are required. From the SINCERA documentation:

    expressions: the FPKM values of 36188 Ensembl genes in 148 cells
    cells: the information of 148 cells, the cluster membership of cells used in the manuscript is encoded in the column "CLUSTER"
    genes: the information of 36188 Ensembl genes
    mouse.ribosomal.genes: a list of ribosomal genes for determining a threshold for specificity filter
    associations.01112014: processed cell type and gene association data downloaded from EBI Expression Atlas (http://www.ebi.ac.uk/gxa/)
    Nkx2.1_data: data for demonstrating the consensus-maximization-based refinement of regulatory target prediction for Nkx2-1 in epithelial cells

This script is for generating these datasets prior to running the pipeline. 
'''

import os,sys,argparse

parser = argparse.ArgumentParser()
parser.add_argument("expressions")
parser.add_argument("cells_file")
parser.add_argument("genes_file")
args = parser.parse_args()

def create_cells(expressions, cells_file):

    clusters = [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 1, 1, 1,
        1, 3, 1, 1, 1, 1, 3, 1, 2, 1, 1, 1, 3, 3, 3, 1, 1, 3, 3, 3, 1, 2, 3, 
        1, 3, 2, 3, 3, 2, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 
        1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
    ]

    clusters_nd = [
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,
        1,1,1,1,2,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,2,2,2,1,2,1,1,2,1,1,
        2,1,1,1,1,2,2,2,1,2,2,2,2,1,2,2,2,1,1,1,2,2,2,2,2,2,2,2,1,2,2,1,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2
    ]

    one_clusts = [1 for x in range(len(clusters))]

    with open(expressions, "r") as ex:
        with open(cells_file, "w") as c:

            for line in ex:
                # print(line)
                cols = line.rstrip("\n").split("\t")
                break
            c.write("CELL\tSAMPLE\tCLUSTER\n")
            count = 0
            for col in cols[1:]:
                # c.write(col + "\t" + col + "\t1\t"+str(clusters[count])+"\n")
                c.write(col + "\t" + col + "\t1\t"+str(clusters_nd[count])+"\n")
                # c.write(col + "\t" + col + "\t1\t"+str(one_clusts[count])+"\n")
                count += 1

def create_genes(expressions, genes_file):

    with open(expressions, "r") as ex:
        with open(genes_file, "w") as g:
            row = 0
            g.write("Geneid\tSYMBOL\tTF\tTF_SOURCE\n")
            for line in ex:
                if row == 0:
                    row += 1
                    continue
                cols = line.rstrip("\n").split("\t")
                gene_entry = cols[0]
                ge_split = gene_entry.split(";")
                gene = ge_split[0]
                g.write(gene_entry + "\t" + gene + "\t1\t1\n")





create_cells(args.expressions, args.cells_file)
create_genes(args.expressions, args.genes_file)