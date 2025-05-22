#Allele-frequency based calculator for 2-3 way genetic admixture calculation 
#23andme and AncestryDna parsed files are all in -- > 



1. Prepare Your Data
23andMe File: Ensure you have a 23andMe file (e.g., example_23andme.txt).

Reference Population Files: Prepare CSV files for each reference population. Each file should have the following columns:

rsid: SNP identifier.

ref_allele: Reference allele.

alt_allele: Alternate allele.

alt_freq: Frequency of the alternate allele in the population.

Example reference population file (pop1.csv):


rsid,ref_allele,alt_allele,alt_freq
rs123456,A,G,0.25
rs234567,C,T,0.10

2. Create a Python Script
Save the AdmixtureCalculator class in a Python script (e.g., admixture_calculator.py). Then, write a script to use the class.

Example usage script (main.py):



from admixture_calculator import AdmixtureCalculator

def main():
    calculator = AdmixtureCalculator()

    calculator.parse_23andme("example_23andme.txt")

    calculator.load_reference_population("Population1", "pop1.csv")
    calculator.load_reference_population("Population2", "pop2.csv")
    calculator.load_reference_population("Population3", "pop3.csv")

    common_snps = calculator.find_common_snps()
    print(f"Number of common SNPs: {len(common_snps)}")

    proportions = calculator.estimate_ancestry()
    print("Ancestry Proportions:")
    for pop, prop in proportions.items():
        print(f"{pop}: {prop:.4f}")

    calculator.visualize_results(proportions)

if __name__ == "__main__":
    main()
3. Run the Script
Save the AdmixtureCalculator class in admixture_calculator.py.

Save the usage script in main.py.

Place your 23andMe file and reference population files in the same directory as the scripts.

Run the script:


python main.py
4. Expected Output
The script will:

Parse the 23andMe file.

Load the reference population data.

Find common SNPs between the 23andMe file and the reference populations.

Estimate ancestry proportions using a likelihood-based approach.

Display the ancestry proportions as a pie chart.

Example output:>>


Number of common SNPs: 500000
Ancestry Proportions:
Population1: 0.4500
Population2: 0.3000
Population3: 0.2500
5. Debugging and Troubleshooting
Missing SNPs: If the number of common SNPs is too low, ensure the reference population files include the same SNPs as the 23andMe file.

File Paths: Double-check the file paths for the 23andMe file and reference population files.

Dependencies: Ensure all required libraries (pandas, numpy, scipy, matplotlib) are installed:

pip install pandas numpy scipy matplotlib
6. Package into a .exe File
If you want to distribute this as a standalone .exe file, use PyInstaller as described earlier. Run the following command:


pyinstaller --onefile main.py
The .exe file will be created in the dist folder.

