import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
import logging
from concurrent.futures import ThreadPoolExecutor
from functools import partial

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AdmixtureCalculator:
    def __init__(self):
        self.reference_pops = {}
        self.individual_snps = None
        self.common_snps = None
        
    def parse_23andme(self, file_path: str) -> pd.DataFrame:
        data = pd.read_csv(file_path, 
                          comment='#', 
                          delim_whitespace=True,
                          names=['rsid', 'chromosome', 'position', 'genotype'])
        data['position'] = pd.to_numeric(data['position'], errors='coerce')
        self.individual_snps = data.dropna()
        return self.individual_snps
    
    def load_reference_population(self, pop_name: str, file_path: str) -> None:
        df = pd.read_csv(file_path)
        if not {'rsid', 'ref_allele', 'alt_allele', 'alt_freq'}.issubset(df.columns):
            raise ValueError("Missing required columns")
        self.reference_pops[pop_name] = df.set_index('rsid')
    
    def find_common_snps(self) -> np.ndarray:
        if self.individual_snps is None or len(self.reference_pops) < 3:
            raise ValueError("Insufficient data")
            
        common_rsids = set(self.individual_snps['rsid'])
        for pop_df in self.reference_pops.values():
            common_rsids &= set(pop_df.index)
            
        self.common_snps = np.array(list(common_rsids))
        return self.common_snps
    
    @staticmethod
    @np.vectorize
    def calculate_genotype_likelihood(genotype: str, ref_freq: float) -> float:
        if len(genotype) != 2:
            return 0.0
        if genotype[0] == genotype[1]:
            return (1 - ref_freq) ** 2 if genotype[0] == 'A' else ref_freq ** 2
        return 2 * ref_freq * (1 - ref_freq)
    
    def estimate_ancestry(self) -> Dict[str, float]:
        def objective(props: np.ndarray) -> float:
            if not np.isclose(np.sum(props), 1.0):
                return np.inf
                
            snp_data = self.individual_snps[self.individual_snps['rsid'].isin(self.common_snps)]
            likelihoods = np.zeros((len(self.common_snps), len(self.reference_pops)))
            
            for i, (pop_name, pop_df) in enumerate(self.reference_pops.items()):
                pop_freqs = pop_df.loc[self.common_snps, 'alt_freq'].values
                genotypes = snp_data.set_index('rsid').loc[self.common_snps, 'genotype'].values
                likelihoods[:, i] = self.calculate_genotype_likelihood(genotypes, pop_freqs)
            
            weighted_likes = np.sum(likelihoods * props, axis=1)
            return -np.sum(np.log(weighted_likes + 1e-300))
        
        n_pops = len(self.reference_pops)
        result = minimize(
            objective,
            x0=np.ones(n_pops) / n_pops,
            bounds=[(0, 1)] * n_pops,
            constraints={'type': 'eq', 'fun': lambda x: np.sum(x) - 1},
            method='SLSQP'
        )
        
        if not result.success:
            raise RuntimeError("Optimization failed to converge")
            
        return dict(zip(self.reference_pops.keys(), result.x))
    
    def visualize_results(self, proportions: Dict[str, float]) -> None:
        plt.figure(figsize=(10, 8))
        plt.pie(
            proportions.values(),
            labels=proportions.keys(),
            autopct='%1.1f%%',
            colors=plt.cm.Set3(np.linspace(0, 1, len(proportions)))
        )
        plt.title('Ancestry Proportions')
        plt.axis('equal')
        plt.show()

