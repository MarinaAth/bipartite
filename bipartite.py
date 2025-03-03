import pandas as pd
from scipy.stats import mannwhitneyu
import os

# Import the binary files
path_binary_df = pd.read_csv('grouped_pathways_binary.tsv',
                             sep='\t', header=0)
path_binary_df.set_index('Pathways', inplace=True)
temp_path = path_binary_df.copy()
path_binary_df.columns = [str(col) + '_R1' for col in path_binary_df.columns]
temp_path.columns = [str(col) + '_R2' for col in temp_path.columns]
merged_path = pd.concat([path_binary_df, temp_path], axis=1)

path_binary_dict = merged_path.transpose().to_dict()

# Import the expression matrix
gene_expression_df = pd.read_csv('normalized_counts_Deseq_Nov2023.data',
                                 sep='\t', header=0)
gene_expression_df.set_index('Gene', inplace=True)
gene_expression_dict = gene_expression_df.to_dict('index')

modDir = ("/Users/mathanasouli/Documents/newRNAseq_Nov2023/"
          "test_bipartite_code")

# Create a dictionary for the interactions that will be saved as s 2d matrix:
# Keys (rows): Modules
# Values (columns): Pathway groups
interaction_dir = {}

for k, v in path_binary_dict.items():
    presence = []
    absence = []
    for k1, v1 in v.items():
        if v[k1] == 1:
            presence.append(k1)
        else:
            absence.append(k1)
    if not presence:
        print(f"Group of pathways {k} is absent in all strains")
    elif not absence:
        print(f"Group of pathways {k} is present in all strains")
    else:
        # Opening the top module files one by one
        print(f"Processing group of pathways {k}")
        for filename in os.listdir(modDir):
            with open(os.path.join(modDir, filename), mode='r') as f:
                # Store the genes in a list
                print(f"Entering {filename}")
                gene_lst = []
                p_temp = []
                a_temp = []
                p_expression = []
                a_expression = []
                gene_sign = []
                interactions = {}
                while line := f.readline():
                    gene_lst.append(line.rstrip())
                    # Calculate the number of genes in the module:
                n = len(gene_lst)
                print(f"Read gene list with length {n}")
                # Init the number of statistically significant genes: N1
                n1 = 0
                # For each gene in the list
                for gene in gene_lst:
                    # Match gene + strains + pattern to gene_expression_df
                    if gene in gene_expression_dict.keys():
                        gene_key = gene_expression_dict[gene]
                        for k_expr, v_expr in gene_key.items():
                            if k_expr in presence:
                                p_expression.append(gene_key[k_expr])
                            elif k_expr in absence:
                                a_expression.append(gene_key[k_expr])
                            else:
                                continue
                            if p_expression and a_expression:
                                statistic, p_value = mannwhitneyu(p_expression, a_expression)
                                if p_value < 0.05:
                                    # +1 counter at n1
                                    n1 += 1
                                    # Add gene to list
                                    gene_sign.append(gene)
                            # Calculate N1/N and assign to interaction
                    else:
                        print(f"Gene {gene} has low expression/excluded.")
                        continue
                mod_path_interact = n1 / n
                # TODO: Check dictionary for more files and save to a file
                interactions[filename] = mod_path_interact
                # Add to dictionary
                interaction_dir[k] = interactions
                # Save gene list to file
                # with open(filename, 'w') as fp:
                #     fp.write("\n".join(str(item) for item in gene_sign))
print(interaction_dir)
