import os, sys
import numpy as np

'''
Authors: Saideep Gona & Martin Ma

This script is for pre-processing single-cell RNA-seq data for the 
course 02715-Advanced Topics in Computational Genomics

The data source comes from a technique called "Mosaic-Seq" which
is a CRISPR-based high-throughput screening metho. It perturbs 
enhancer function in order to allow for downstream analysis of 
expression changes.
'''

# Assuming input is an expression table of cells x genes

class SingleCellAnalysis():

    def __init__(self, 
                filename, 
                filter_criteria, 
                imputation_method, 
                normalization_method,
                clustering_method,
                diff_genes_method):
        
        self.load_data(filename)

        self.steps

        self.filter_method = filter_method 
        self.imputation_method = imputation_method
        self.normalization_method = normalization_method
        self.clustering_method = clustering_method
        self.diff_genes_method = diff_genes_method

        self.pp_steps = [
            [self.filter_table, self.filter_method],
            [self.impute_dropout, self.imputation_method],
            [self.normalize, self.normalization_method],
            [self.cluster_cells, self.clustering_method],
            [self.diff_genes, self.diff_genes_method]
        ]

    def run_preprocessing(self):
        '''
        Runs suite of preprocessing steps in sequence
        '''
        
        for step in self.pp_steps:


    def load_data(self, filename):
        '''
        Load the source data
        '''
        loaded = np.loadtxt(filename, ",")
        self.start_table = loaded

    def filter_table(self):
        '''
        Function for filtering out unwanted data 
        '''

        self.imputed_table = imputed_table

    def impute_dropout(self):
        '''
        Function for imputing dropout given an expression matrix
        '''

        self.imputed_table = imputed_table

    def normalize(self):
        '''
        Function for normalizing expression data prior to further 
        analysis
        '''

        self.normalized_table = normalized_table

    def cluster_cells(self):
        '''
        Perform clustering on the expression matrix. 
            Input: cell x gene expression table
            Output: Clusters of cells
        '''

        self.clusters = clusters

    def diff_genes(self):
        '''
        Computes differentially expressed genes between a set of clusters
        based on the expression table
        '''

        self.diff_genes = diff_genes

input_data_file = "mosaic_seq_processed.csv"

start_table = load_data(input_data_file)

filtered_table = filter_table(start_table)

imputed_table = impute_dropout(filtered_table)

normalizaed_table = 