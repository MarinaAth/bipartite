import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as multitest
import save_mods_to_dict
import numpy as np


modDir = ('/Users/mathanasouli/Documents/final_network_Sept2024/topMods')
modules_all = save_mods_to_dict.assign_to_dict(modDir)
process_mods = {f'mod{i}' for i in range(50, 61)}
modules = {k: modules_all[k] for k in process_mods}
print(modules.keys(), flush=True)

# Import the binary files
path_binary_df = pd.read_csv('pathways_to_test',
                             sep='\t', header=0)
path_binary_df.set_index('Pathways', inplace=True)
temp_path = path_binary_df.copy()
path_binary_df.columns = [str(col) + '_R1' for col in path_binary_df.columns]
temp_path.columns = [str(col) + '_R2' for col in temp_path.columns]
merged_path = pd.concat([path_binary_df, temp_path], axis=1)
merged_path.to_csv("grouped_pathways_reps.tsv", sep="\t")

path_binary_dict = merged_path.transpose().to_dict()

# Import the expression matrix
gene_expression_df = pd.read_csv('normalized_counts_final_Sep2024.tsv',
                                 sep='\t', header=0)
gene_expression_df.set_index('Gene', inplace=True)
gene_expression_dict = gene_expression_df.to_dict('index')

# Create a dictionary for the interactions that will be saved as s 2d matrix:
# Keys (rows): Modules
# Values (columns): Pathway groups
interaction_dir_FDR = {}
interaction_dir_Bonferroni = {}

p_temp = []
a_temp = []
presence_expr = []
absence_expr = []
gene_lst = []
tmp_FDR = {}
tmp_Bonferroni = {}
mod_path_interact_FDR = 0
mod_path_interact_Bonferroni = 0
gene_pvals = []
sig_genes = []
for k, v in path_binary_dict.items():
    presence = []
    absence = []
    for k1, v1 in v.items():
        if v[k1] == 1:
            presence.append(k1)
        else:
            absence.append(k1)
    if not presence:
        print(f"Group of pathways {k} is absent in all strains", flush=True)
        continue
    elif not absence:
        print(f"Group of pathways {k} is present in all strains", flush=True)
        continue
    else:
        print(f"Processing group of pathways {k}", flush=True)
        # Opening the top module files one by one
        for k_mod, v_mod in modules.items():
            gene_lst = modules[k_mod]
            # Calculate the number of genes in the module:
            n = len(gene_lst)
            # Init the statistically significant genes:
            # N1 for FDR and Bonferroni
            n1_FDR = 0
            n1_Bonferroni = 0
            # For each gene in the list
            print(f"Evaluating  {k_mod} against {k}", flush=True)
            gene_pvals.clear()
            sig_genes.clear()
            p_temp.clear()
            a_temp.clear()
            presence_expr.clear()
            absence_expr.clear()
            for gene in gene_lst:
                # Match gene + strains + pattern to gene expression
                if gene in gene_expression_dict.keys():
                    gene_key = gene_expression_dict[gene]
                    gene_expr_values = []
                    for sublist in gene_key.values():
                        if isinstance(sublist, list):
                            # If it's a list, extend
                            gene_expr_values.extend(sublist)
                        else:
                            gene_expr_values.append(sublist)
                    median_gene = np.median(gene_expr_values)
                    for k_expr, v_expr in gene_key.items():
                        if k_expr in presence:
                            presence_expr.append(gene_key[k_expr])
                        elif k_expr in absence:
                            absence_expr.append(gene_key[k_expr])
                        else:
                            continue
                    # Calculate p-value and add to list
                    fc_present = [element / median_gene for element in presence_expr]
                    fc_absent = [element / median_gene for element in absence_expr]
                    fc_present_median = np.median(fc_present)
                    fc_absent_median = np.median(fc_absent)
                    if (fc_present_median > (fc_absent_median * 2) and
                            presence_expr and absence_expr):
                        statistic, p_value = stats.mannwhitneyu(presence_expr,
                                                                absence_expr)
                        gene_pvals.append(p_value)
                    else:
                        tmp_FDR[k_mod] = 0.0
                        tmp_Bonferroni[k_mod] = 0.0
                else:
                    print(f"Gene {gene} has low expression/excluded.",
                          flush=True)
                    continue
            # Sort list of p-values
            gene_pvals.sort()
            # Perform FDR and Bonferroni correction
            if gene_pvals:
                gene_pvals.sort()
                bonferroni_test = np.array(gene_pvals) * len(gene_pvals)
                reject, fdr_test, _, _ = multitest.multipletests(gene_pvals,
                                                           alpha=0.05,
                                                           method='fdr_bh')
                # Filter genes with corrected values less than 0.05
                sig_FDR = [x for x in fdr_test if x < 0.05]
                sig_Bonferroni = [y for y in bonferroni_test if y < 0.05]
                # Calculate interacting genes in both correction methods
                n1_FDR = len(sig_FDR)
                n1_Bonferroni = len(sig_Bonferroni)
                # Calculate interaction values
                mod_path_interact_FDR = (n1_FDR / n)
                mod_path_interact_Bonferroni = (n1_Bonferroni / n)
                # Update the temporary dictionaries
                tmp_FDR[k_mod] = float(mod_path_interact_FDR)
                tmp_Bonferroni[k_mod] = float(mod_path_interact_Bonferroni)
            else:
                continue
    # Add to respective dictionaries
    interaction_dir_FDR[k] = tmp_FDR
    interaction_dir_Bonferroni[k] = tmp_Bonferroni

# Write interaction with FDR to file
bipartite_df_FDR = pd.DataFrame.from_dict(interaction_dir_FDR)
bipartite_df_FDR.to_csv("bipartite_mod50_60_FDR_Sep2024.tsv", sep="\t")

# Write interaction with Bonferroni to file
bipartite_df_Bonferroni = pd.DataFrame.from_dict(interaction_dir_Bonferroni)
bipartite_df_Bonferroni.to_csv("bipartite_mod50_60_Bonferroni_Sep2024.tsv", sep="\t")