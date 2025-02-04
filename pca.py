import pandas as pd
#simple python script to follow through 

def load_23andme_data(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    df.columns = ['rsid', 'chromosome', 'position', 'genotype']
    
    df = df[df['genotype'].str.match(r'^[ATCG]{2}$')]
    
    df['genotype'] = df['genotype'].apply(lambda x: sum([{'A': 0, 'T': 1, 'C': 2, 'G': 3}[base] for base in x]))
    
    return df

def load_reference_data(file_path):
    df = pd.read_csv(file_path)
    
    if 'rsid' not in df.columns or 'genotype' not in df.columns:
        raise ValueError("Reference dataset must contain 'rsid' and 'genotype' columns.")
    
    return df

def combine_data(data_23andme, reference_data):
    combined_data = pd.merge(data_23andme, reference_data, on='rsid', suffixes=('_23andme', '_ref'))
    
    return combined_data

def save_combined_data(combined_data, output_file):
    combined_data.to_csv(output_file, index=False)

file_path = 'your_23andme_file.txt'
reference_file_path = 'reference_dataset.csv'
output_file = 'combined_data.csv'

data_23andme = load_23andme_data(file_path)
reference_data = load_reference_data(reference_file_path)

combined_data = combine_data(data_23andme, reference_data)

save_combined_data(combined_data, output_file)

print(f"Combined data saved to {output_file}")
