## SYNOPSIS


PACtool is a programme that performs GWAS analyses on control/case data of SNP variants and was implemented using Python 3.5

PACtool 1.x Features:
   
1. Allele frequency calculation of each variant for all selected SNPs
2. HWE statistic and  respective p-value calculation for the unified controls and cases dataset
3. Linkage Disequilibrium evaluation between two selected SNPs  (D' and r-squared are calculated)
4. Association Test on the provided dataset, with optional generation of manhattan plot and qq-plot
5. Chromosome coordinates conversion from ncbi 36 assembly (hg18) to 38 (hg38)
6. Information retrieval for the selected SNPs, from Ensembl's Variant Effect Predictor database.

PACtool Limitations:
	 
1. Current version of PACtool does not support analysis for variants located on different chromosomes.   
   Please make sure that all SNPs in your dataset are located on the same chromosome.


## INPUT FILES FORMAT

The input files must be provided in accordance to the HAPGEN2 programme, which can be
applied on various genomic datasets. The input file consists of rows and columns, 
where each row represents a SNP, and the columns contain information about the following: 

Column_1: snp_id
	-e.g. snp_0
	-a unique identifier for each SNP 
	

Column_2: rs_id 
	-e.g. rs6054257
	-or alternatively  the genomic coordinates, with a prefix indicating the chromosome e.g. 20-9150

Column_3: snp_coordinates based on NCBI build 36
	-e.g. 9150

Column_4: reference allele
	-denoted with the first letter of each nucleotide, e.g. A if the reference allele has Adenine in the point coordinates indicated in column 3

Column_5: alternative allele
	-denoted with the first letter of each nucleotide as well, e.g. G if the alternate allele has Guanine in the point coordinates indicated in column 3

Column_6 – Last_Column:
	-These columns contain the genotype for each sample. 
	 Each sample genotype is represented by 3 columns consisted of '0' and '1' digits, 
	 with the positioning of '1' being the indicator of the genotype. 
	 More specifically,
   '''
	- 1 0 0 --->  ref-ref , homozygous for the reference allele 
	- 0 1 0 --->  ref-alt , heterozygous
	- 0 0 1 --->  alt-alt , homozygous for the alternative allele
  '''

## OUTPUT FILES FORMAT

The output files follow the same white space separated file format as the input files described above. 
The first column for each SNP in all output files is the snp_id as described above, 
while the rest of the columns hold the values from the respective statistical tests and analyses.


## INSTALLATION-UTILIZATION

PACtool was uploaded to the Python Package Index (PyPI), using the standard procedure 
(creating source distribution, registering the package against PyPI server, and uploading it). 

PACtool can be downloaded by visiting the following link:

	https://pypi.python.org/pypi/pactool

The package will be downloaded as a .tgz file, which you can subsequently untar and unzip. 
In the pactool directory, a message with the necessary information on how to execute the programme will appear by typing:

	python3 pactool.py -h

You can execute the file pactool.py using python 3, and include all preferred arguments. 
Three of the arguments are required and must be provided every time, or else you will receive an error message. 
Please be sure to include the following arguments:

	-controls_file: Indicates the input file containing control samples
	-cases_file:	Indicates the input file o
	-output:	Specifies the prefix of each output file

The following optional arguments can be also include, alongside with a file with snp_codes (one in each line) to 
perform the corresponding action:
	
	-keep_snps:	Keeps for analysis only the SNPs specified in the provided file.
	-remove_snps:	Removes from following analysis the SNPs specified in the provided file.

The above actions can be applied to the samples as well (the lines in given file shall be e.g. control_5 or case_10):

	-keep_samples:   Keeps for analysis only the control/case samples specified in the provided file.
	-remove_samples: Removes from further analyisis the control/case samples specified in the provided file.

After performing the keep/remove actions, the following analysis options can be selected:

	-allele_frequency 	Calculates the frequencies of the reference and alternative variants in control samples, case samples as well as their total frequencies.
				Outputs the file 'output'.frequency with 7 columns: snp_code ref_freq_control alt_freq_control ref_freq_cases alt_freq_cases ref_freq_total alt_freq_total
	-hwe			Calculates the Hardy-Weinberg Equilibrium statistic and the corresponding p-value.
				Outputs the file 'output'.hwe with 3 columns: snp_code hwe_statistics p-value
	-ld SNP1 SNP2		Estimates if the two given SNPs are in Linkage Disequilbrium, by calculating the D' and r-squared statistics. SNP1 and SNP2 are required snp_codes.
				Outputs the file 'output'.ld with 4 columns: snp1_code snp2_code D' r-squared
	-association_test	Performs genotypic association test for each SNP and calculates the odds-ratios (r=reference, a=alternative)
				Outputs the file 'output'.association with 8 columns: snp_code locus ref alt p-value OR_rr_ra OR_rr_aa OR_ra_aa
	-manhattan		Draws a manhattan plot for the p-values of the association test.
				Can only be used if -association_test argument is given.
	-qqplot			Draws a qq-plot for the p-values of the association test.
				Can only be used if -association_test argument is given.
	-get_info SNP		Retrieves information about variant with snp_code SNP. Prints a json format output with all
				information obtained from VEP database of Ensembl.



## TESTING

PACtool was built around and tested using the following files:

gwas.controls.gen
gwas.cases.gen


which can be download from the following link:

https://s3.eu-central-1.amazonaws.com/pythonprojectgwas/gwas.tar.gz

Credits to the creation of this artificial dataset goes to the team running 
'BIO-102 Introduction to Programming' course of MSc Bioinformatics, UoCrete, Heraklion.



## CONTRIBUTORS

Authors of this project are:

Christina Chatzipantsiou  (chatzipantsiou@gmail.com)
Panayiotis Linardos       (mondestrasz@gmail.com)
Paschalis Natsidis        (pnatsidis@hotmail.com)

For any bug report, contribution or general comment please contact any of the authors on the provided e-mail addreses.


## LICENSE

Pactool lies under MIT license described in this link:

https://opensource.org/licenses/MIT
