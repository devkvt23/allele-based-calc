import pandas as pd

def load_23andme_data(file_path):
    # Load the 23andMe data
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    df.columns = ['rsid', 'chromosome', 'position', 'genotype']
    
    # Filter out non-SNP data (e.g., indels, no-calls)
    df = df[df['genotype'].str.match(r'^[ATCG]{2}$')]
    
    # Convert genotypes to numerical values (e.g., AA=0, AG=1, GG=2)
    df['genotype'] = df['genotype'].apply(lambda x: sum([{'A': 0, 'T': 1, 'C': 2, 'G': 3}[base] for base in x]))
    
    return df

def load_reference_data(file_path):
    # Load the reference dataset
    df = pd.read_csv(file_path)
    
    # Ensure the reference dataset has the required columns
    if 'rsid' not in df.columns or 'genotype' not in df.columns:
        raise ValueError("Reference dataset must contain 'rsid' and 'genotype' columns.")
    
    return df

def combine_data(data_23andme, reference_data):
    # Merge the datasets on 'rsid'
    combined_data = pd.merge(data_23andme, reference_data, on='rsid', suffixes=('_23andme', '_ref'))
    
    return combined_data

def save_combined_data(combined_data, output_file):
    # Save the combined data to a file
    combined_data.to_csv(output_file, index=False)

# Example usage
file_path = 'your_23andme_file.txt'
reference_file_path = 'reference_dataset.csv'
output_file = 'combined_data.csv'

# Load data
data_23andme = load_23andme_data(file_path)
reference_data = load_reference_data(reference_file_path)

# Combine data
combined_data = combine_data(data_23andme, reference_data)

# Save combined data
save_combined_data(combined_data, output_file)

print(f"Combined data saved to {output_file}")