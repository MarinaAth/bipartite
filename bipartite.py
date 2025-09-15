import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as multitest
import save_mods_to_dict
import numpy as np
from typing import Dict, List, Tuple, Optional
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class BipatiteAnalysis:
    """Handles the construction of a bipartite network between
    co-expression modules (nematode) and present/absent bacterial pathways
    """

    def __init__(self, mod_dir:str, pathway_file:str, expression_file:str,
                 module_range: Tuple[int, int] = (50,61)):
        self.mod_dir = mod_dir
        self.pathway_file = pathway_file
        self.expression_file = expression_file
        self.module_range = module_range

        # Initialising input structures
        self.modules = {}
        self.path_binary_dict = {}
        self.gene_expression = {}

        # Initialising result structures
        self.interaction_results_fdr = {}
        self.interaction_results_bonferroni = {}


    def load_modules(self) -> None:
        """Importing and processing bacterial pathway
        data (binary for presence/absence)
        """
        logger.info(f"Loading pathways from {self.pathway_file}")
        path_binary_df = pd.read_csv(self.pathway_file, sep='\t', header=0)
        path_binary_df.set_index('Pathways', inplace=True)

        # Duplicate columns for each of 2 replicates
        temp_path = path_binary_df.copy()
        path_binary_df.columns = [str(col) + '_R1' for col in path_binary_df.columns]
        temp_path.columns = [str(col) + '_R2' for col in temp_path.columns]

        # Merge and save
        merged_path = pd.concat([path_binary_df, temp_path], axis=1)
        output_file = "grouped_pathways_reps.tsv"
        merged_path.to_csv(output_file, sep="\t")
        logger.info(f"Created and merged binary pathway file for replicates")
        self.path_binary_dict = merged_path.transpose().to_dict()


    def load_gene_expression(self) -> None:
        """Importing and processing the gene expression data
        from the nematode's transcriptomic profiles
        """
        logger.info(f"Loading gene expression table")
        gene_expression_df = pd.read_csv(self.expression_file, sep='\t', header=0)
        gene_expression_df.set_index('Gene', inplace=True)
        self.gene_expression_dict = gene_expression_df.to_dict('index')
        logger.info(f"Loaded gene expression data for {len(self.gene_expression_dict)} genes")


    def get_presence_absence_strains(self, pathway_data: Dict) -> Tuple[List[str], List[str]]:
        """Separate strains based on the presence/absence patterns in the pathway data

        Args:
            pathway_data (Dict): 

        Returns:
            Tuple[List[str], List[str]]: _description_
        """
        presence = [strain for strain, value in pathway_data.items() if value == 1]
        absence = [strain for strain, value in pathway_data.items() if value == 0]
        return presence, absence
    

    def calculate_fold_change_pval(self, gene:str, presence_strains: List[str],
                                   absence_strains: List[str]) -> Optional[float]:
        """Calculate fold change and pvalue of FC for each gene of each module
        for strains with present vs. absent pathways

        Args:
            gene (str): _description_
            presence_strains (List[str]): _description_
            absence_strains (List[str]): _description_

        Returns:
            Optional[float]: _description_
        """
        if gene not in self.gene_expression_dict:
            return None
        
        gene_data = self.gene_expression_dict[gene]

        # Extract expression values of gene and ensure there are expression values 
        # for both presence and absence (for FC and pval)
        presence_expr = [gene_data[strain] for strain in presence_strains if strain in gene_data]
        absence_expr = [gene_data[strain] for strain in absence_strains if strain in gene_data]

        if not presence_expr or not absence_expr:
            return None
        
        # Median expression across all diets
        all_expr_values = list(gene_data.values())
        median_gene = np.median(all_expr_values)

        # FC calculation
        fc_present = [expr / median_gene for element in presence_expr]
        fc_absent = [expr / median_gene for element in absence_expr]

        fc_present_median = np.median(fc_present)
        fc_absent_median = np.median(fc_absent)

        # Test only for FC threshold
        if fc_present_median > (fc_absent_median * 2):
            try:
                statistic, p_value = stats.mannwhitneyu(presence_expr, absence_expr)
                return p_value
            except ValueError as e:
                logger.warning(f"Mann-Whitney U test failed for gene {gene}: {e}")
                return None
            

    def analyze_module_pathway_interaction(self, module_genes: List[str],
                                           presence_strains: List[str],
                                           absence_strains: List[str]) -> Tuple[float, float]:
        """Interaction between a gene co-expression module and a bacterial pathway group

        Args:
            module_genes (List[str]): _description_
            presence_strains (List[str]): _description_
            absence_strains (List[str]): _description_

        Returns:
            Tuple[float, float]: _description_
        """
        gene_pvals = []

        for gene in module_genes:
            p_value = self.calculate_fold_change_pval(gene, presence_strains, absence_strains)
            if p_value is not None:
                gene_pvals.append(p_value)

        if not gene_pvals:
            return 0.0, 0.0
        
        # Multiple testing correction
        gene_pvals.sort()
        n_genes = len(module_genes)

        # Bonferroni
        bonferroni_corrected = np.array(gene_pvals) * len(gene_pvals)
        sig_Bonferroni = np.sum(bonferroni_corrected < 0.05)

        # FDR
        reject, fdr_corrected, _, _ = multitest.multipletests(gene_pvals, alpha=0.05, method='fdr_bh')
        sig_fdr = np.sum(fdr_corrected < 0.05)

        # Calculate interaction scores
        fdr_score = sig_fdr / n_genes
        bonferroni_score = sig_Bonferroni / n_genes

        return fdr_score, bonferroni_score
    

    def run_analysis(self) -> None:
        """Run bipartite
        """
        logger.info("Starting bipartite")

        self.load_modules()
        self.load_pathways()
        self.load_gene_expression()

        # Process each bacterial MPG
        for pathway_name, pathway_data in self.path_binary_dict.items():
            presence_strains, absence_strains = self.get_presence_absence_strains(pathway_data)

            if not presence_strains:
                logger.info(f"MPG {pathway_name} is absent in all strains")
                continue
            elif not absence_strains:
                logger.info(f"MPG {pathway_name} is present in all strains")
                continue
            
            logger.info(f"Processing MPG: {pathway_name}")

            # Initialise results for MPG
            pathway_results_fdr = {}
            pathway_results_bonferroni = {}

            # Process each module
            for module_name, module_genes in self.modules.items():
                logger.info(f"Evaluating {module_name} against {pathway_name}")

                fdr_score, bonferroni_score = self.analyze_module_pathway_interaction(module_genes,
                                                                                      presence_strains,
                                                                                      absence_strains)
                pathway_results_fdr[module_name] = fdr_score
                pathway_results_bonferroni[module_name] = bonferroni_score

            self.interaction_results_fdr[pathway_name] = pathway_results_fdr
            self.interaction_results_bonferroni[pathway_name] = pathway_results_bonferroni


    def save_results(self, output_prefix:str = 'bipartite_mod50_60') -> None:
        fdr_df = pd.DataFrame.from_dict(self.interaction_results_fdr)
        fdr_output = f"{output_prefix}_FDR.tsv"
        fdr_df.to_csv(fdr_output, sep="\t")
        logger.info(f"FDR results saved to {fdr_output}")

        bonferroni_df = pd.DataFrame.from_dict(self.interaction_results_bonferroni)
        bonferroni_output = f"{output_prefix}_Bonferroni.tsv"
        bonferroni_df.to_csv(bonferroni_output, sep="\t")
        logger.info(f"Bonferroni results saved to {bonferroni_output}")


def main():
    config = {
        'mod_dir': 'final_network/topMods',
        'pathway_file': 'binary_pathways',
        'expression_file': 'normalized_counts.tsv',
        'module_range': (50, 61),
        'output_prefix': 'bipartite_mod_50_61',
    }

    analyzer = BipatiteAnalysis(
        mod_dir=config['mod_dir'],
        pathway_file=config['pathway_file'],
        expression_file=config['expression_file'],
        module_range=config['module_range']
    )

    analyzer.run_analysis()
    analyzer.save_results(config['output_prefix'])

    logger.info("Completed..")


if __name__ == "__main__":
    main()

