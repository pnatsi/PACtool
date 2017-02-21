#!/usr/bin/python3

import argparse
import frequencies
import kor
import timeit
import hwe
import getinfo
import ld
import association_test
import multitools

start = timeit.default_timer()

#HERE STARTS THE ARGUMENT DEFINING
usage = "A program for GWAS studies..."
toolname = "pactool"
footer = "Who \n Christina Chatzipantsiou (chatzipantsiou@gmail.com);  \n \
Panayiotis Linardos (mondestrasz@gmail.com); \n \
Paschalis Natsidis (pnatsidis@hotmail.com); \n \nWhere \n Introduction to Programming;\n\
 MSc Bioinformatics, University of Crete; \n  \nWhen\n February 2017; \n\n"

parser = argparse.ArgumentParser(description = usage, prog = toolname, epilog = footer, formatter_class=argparse.RawDescriptionHelpFormatter,)
parser.add_argument('-controls_file', metavar = 'filename', dest = 'controls', required = True,
                    help = 'file w/ control samples')
parser.add_argument('-cases_file', metavar = 'filename', dest = 'cases', required = True,
                    help = 'file w/ cases samples')
parser.add_argument('-output',  metavar = 'filename', dest = 'output', required = True,
                    help = "name of output file(s), e.g. 'output'.hwe")

group1 = parser.add_mutually_exclusive_group()
group1.add_argument('-keep_snps', metavar = 'filename', dest = 'keep_snps',
                    help = 'file w/ SNP codes. Analyze these SNPs only')
group1.add_argument('-remove_snps',  metavar = 'filename', dest = 'remove_snps',
                    help = 'file w/ SNP codes. Remove these SNPs from analysis')

group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-keep_samples',  metavar = 'filename', dest = 'keep_samples',
                    help = "file w/ sample codes (e.g 'case_5', 'control_100'). Analyze these samples only")
group2.add_argument('-remove_samples',  metavar = 'filename', dest = 'remove_samples',
                    help = "file w/ sample codes (e.g 'case_5', 'control_100'). Remove these samples from analysis")

parser.add_argument('-allele_frequency', action = 'store_true', dest = 'freqs',
                    help = "creates 'output'.frequencies w/ reference and alternative frequencies in control and cases samples")
parser.add_argument('-HWE', '-hwe', action = 'store_true', dest = 'hwe',
                    help = "creates 'output'.hwe w/ Hardy-Weinberg Equilibrium for control+cases dataset")
parser.add_argument('-LD', '-ld', nargs = 2, dest = 'ld',  metavar = ('SNPcode1', 'SNPcode2'),
                    help = "creates 'output'.ld w/ Linkage Disequilibrium (r-squared, D') between two SNPs")
parser.add_argument('-association_test',  action = 'store_true', dest = 'assoc',
                   help = "creates 'output'.association w/ association test results for all SNPs")
parser.add_argument('-manhattan', action = 'store_true', dest = 'manhattan',
                    help = 'draws a manhattan plot for p-values of association test')
parser.add_argument('-qqplot', '-qq', action = 'store_true', dest = 'qqplot',
                    help = 'draws a qq-plot for p-values of association test')
parser.add_argument('-get_info', dest = 'info',  metavar = 'SNPcode', 
                    help = 'pulls information from online databases for this SNP')

#parser.print_help()

args = parser.parse_args()


#READ FILENAMES FROM USER INPUT
controls_file = args.controls
cases_file = args.cases
output_filename = args.output
keepsnps_file = args.keep_snps
removesnps_file = args.remove_snps
keepsamples_file = args.keep_samples
removesamples_file = args.remove_samples

#----KEEP OR REMOVE SNP AND/OR SAMPLES-----------------------------------------------

#if no kor action is defined

controls_to_analyze = []
cases_to_analyze = []

if not keepsnps_file and not removesnps_file and not keepsamples_file and not removesamples_file:
    print("Reading the input files...", end="")
    with open(controls_file, 'r') as clf:
        for line in clf.readlines():
            splitted = line.split()
            controls_to_analyze.append(splitted)
        
    with open(cases_file, 'r') as csf:
        for line in csf.readlines():
            splitted = line.split()
            cases_to_analyze.append(splitted)
    print(" OK!")
    
#keep_snps only
if keepsnps_file and not keepsamples_file and not removesamples_file:
    print("Reading the input files...", end="")
    controls_kor_snps = kor.keep_snps(controls_file, keepsnps_file)
    cases_kor_snps = kor.keep_snps(cases_file, keepsnps_file)
    for i in controls_kor_snps:
        splitted = i.split()
        controls_to_analyze.append(splitted)
    for j in cases_kor_snps:
        splitted = j.split()
        cases_to_analyze.append(splitted)
    print(" OK!")
        
#remove_snps only
if removesnps_file and not keepsamples_file and not removesamples_file:
    print("Reading the input files...", end="")
    controls_kor_snps = kor.remove_snps(controls_file, removesnps_file)
    cases_kor_snps = kor.remove_snps(cases_file, removesnps_file)
    for i in controls_kor_snps:
        splitted = i.split()
        controls_to_analyze.append(splitted)
    for j in cases_kor_snps:
        splitted = j.split()
        cases_to_analyze.append(splitted)
    print(" OK!")
    
#keep_samples only
if keepsamples_file and not keepsnps_file and not removesnps_file:
    print("Reading the input files...", end="")
    kor_samples = kor.keep_samples(controls_file, cases_file, keepsamples_file)
    controls_kor_samples = kor_samples[0]
    cases_kor_samples = kor_samples[1]
    controls_to_analyze = controls_kor_samples
    cases_to_analyze = cases_kor_samples
    print(" OK!")
    
#remove_samples only
if removesamples_file and not keepsnps_file and not removesnps_file:
    print("Reading the input files...", end="")
    kor_samples = kor.remove_samples(controls_file, cases_file, removesamples_file)
    controls_kor_samples = kor_samples[0]
    cases_kor_samples = kor_samples[1]
    controls_to_analyze = controls_kor_samples
    cases_to_analyze = cases_kor_samples
    print(" OK!")

#UNIFIED KOR ACTIONS

#keep_snps and keep_samples
if keepsnps_file and keepsamples_file:
    print("Reading the input files...", end="")
    controls_kor_snps = kor.keep_snps(controls_file, keepsnps_file)
    cases_kor_snps = kor.keep_snps(cases_file, keepsnps_file)
    kor_samples = kor.keep_samples(controls_file, cases_file, keepsamples_file)
    controls_kor_samples = kor_samples[0]
    cases_kor_samples = kor_samples[1]
    for i in controls_kor_snps:
        splitted = i.split()
        for j in controls_kor_samples:
            if j[0] == splitted[0]:
                controls_to_analyze.append(j)
    for i in cases_kor_snps:
        splitted = i.split()
        for j in cases_kor_samples:
            if j[0] == splitted[0]:
                cases_to_analyze.append(j)
    print(" OK!")
                
#remove_snps and remove_samples                
if removesnps_file and removesamples_file:
    print("Reading the input files...", end="")
    controls_kor_snps = kor.remove_snps(controls_file, removesnps_file)
    cases_kor_snps = kor.remove_snps(cases_file, removesnps_file)
    kor_samples = kor.remove_samples(controls_file, cases_file, removesamples_file)
    controls_kor_samples = kor_samples[0]
    cases_kor_samples = kor_samples[1]
    for i in controls_kor_snps:
        splitted = i.split()
        for j in controls_kor_samples:
            if j[0] == splitted[0]:
                controls_to_analyze.append(j)
    for i in cases_kor_snps:
        splitted = i.split()
        for j in cases_kor_samples:
            if j[0] == splitted[0]:
                cases_to_analyze.append(j)
    print(" OK!")     
           
#keep_snps and remove_samples
if keepsnps_file and removesamples_file:
    print("Reading the input files...", end="")
    controls_kor_snps = kor.keep_snps(controls_file, keepsnps_file)
    cases_kor_snps = kor.keep_snps(cases_file, keepsnps_file)
    kor_samples = kor.remove_samples(controls_file, cases_file, removesamples_file)
    controls_kor_samples = kor_samples[0]
    cases_kor_samples = kor_samples[1]
    for i in controls_kor_snps:
        splitted = i.split()
        for j in controls_kor_samples:
            if j[0] == splitted[0]:
                controls_to_analyze.append(j)
    for i in cases_kor_snps:
        splitted = i.split()
        for j in cases_kor_samples:
            if j[0] == splitted[0]:
                cases_to_analyze.append(j)
    print(" OK!")
    
#remove_snps and keep_samples   
if removesnps_file and keepsamples_file:
    print("Reading the input files...", end="")
    controls_kor_snps = kor.remove_snps(controls_file, removesnps_file)
    cases_kor_snps = kor.remove_snps(cases_file, removesnps_file)
    kor_samples = kor.keep_samples(controls_file, cases_file, keepsamples_file)
    controls_kor_samples = kor_samples[0]
    cases_kor_samples = kor_samples[1]
    for i in controls_kor_snps:
        splitted = i.split()
        for j in controls_kor_samples:
            if j[0] == splitted[0]:
                controls_to_analyze.append(j)
    for i in cases_kor_snps:
        splitted = i.split()
        for j in cases_kor_samples:
            if j[0] == splitted[0]:
                cases_to_analyze.append(j)
    print(" OK!")               

#--------------------------------------------------------------------------------------    

if args.freqs:
    print("Calculating frequencies...", end="")
    controls_freq = []
    cases_freq = []
    
    for snp in controls_to_analyze:
        controls_freq.append(frequencies.frequencies(snp))
    for snp in cases_to_analyze:
        cases_freq.append(frequencies.frequencies(snp))
    
    #APPEND CASES FREQUENCIES TO CONTROL FREQUENCIES
    #NOW WE HAVE A LIST WITH 5 ELEMENTS (CODE, REF_CON, ALT_CON, REF_CAS, ALT_CAS)    
    result = []

    for k in range(0, len(controls_freq)):
        current = controls_freq[k]
        current.append(cases_freq[k][1])
        current.append(cases_freq[k][2])
        result.append(controls_freq[k])

    #APPEND TOTAL FREQUENCIES AS 6TH AND 7TH ELEMENT, AS THE AVERAGE OF THE OTHER TWO
    for r in result:
        ref_freq_total = round((r[1]+r[3])/2, 4)
        alt_freq_total = round((r[2]+r[4])/2, 4)
        r.append(ref_freq_total)
        r.append(alt_freq_total)

    #WRITE OUTPUT TO FILE IN 7 COLUMNS SEPARATED WITH SPACE    
    freqs_filename = output_filename + '.frequency'
    with open(freqs_filename, 'w') as output:
        output.write("snp_code" + " " + "ref_freq_control" + " " + "alt_freq_control" + " " + "ref_freq_cases" \
                 + " " + "alt_freq_cases" + " " + "ref_freq_total" + " " + "alt_freq_total" + "\n")
        for r in result:   
            output.write(str(r[0]) + " " + str(r[1]) + " " + str(r[2]) + " " + str(r[3]) + " " +\
                         str(r[4]) + " " + str(r[5]) + " " + str(r[6]) + " " + "\n")
    print(" OK!")  
      
if args.hwe:    
    print("Calculating HWE...", end="")
    united = hwe.hwe(controls_to_analyze, cases_to_analyze)
    hwe_filename = output_filename + '.hwe'
    with open(hwe_filename, 'w') as output:
        output.write("snp_code" + " " + "hwe_statistic" + " " + "p-value" + "\n")
        for entry in united:
            output.write(str(entry[0]) + " " + str(entry[1]) + " " + str(entry[2]) + "\n")
    print(" OK!")
            
if args.ld:
    print("Calculating Linkage Disequilibrium...", end="")
    snp_code1 = args.ld[0]
    snp_code2 = args.ld[1]
    result = ld.ld(snp_code1, snp_code2, controls_to_analyze, cases_to_analyze)
    ld_filename = output_filename + '.ld'
    with open(ld_filename, 'w') as output:
        output.write("snp1_code" + " " + "snp2_code" + " " + "D'" + " " + "r-squared" + "\n")
        output.write(snp_code1 + " " + snp_code2 + " " + str(result[0]) + " " + str(result[1]) + "\n")
    print(" OK!")
    
if args.assoc:
    print("Performing association test...", end="")
    association_filename = output_filename + '.association'
    with open(association_filename, 'w') as output:
            output.write("snp_code" + " " + "locus" + " " + "ref" + " " + "alt" + " " + "p-value" + " " + "OR_rr_ra" + " " + "OR_rr_aa" + " " + "OR_ra_aa" + "\n")
            for i in range(len(controls_to_analyze)): #controls and cases have the same length (even if we remove snps, they are removed from both groups)
                pval, OR = association_test.genotype_test(controls_to_analyze[i], cases_to_analyze[i])
                current = controls_to_analyze[i]
                snpNo=current[0]
                locus=current[2]
                ref=current[3] 
                alt=current[4]
                OR_rr_ra = round(OR[ref+ref+"/"+ref+alt], 4)
                OR_rr_aa = round(OR[ref+ref+"/"+alt+alt], 4)
                OR_ra_aa = round(OR[ref+alt+"/"+alt+alt], 4)
                output.write(snpNo + " " + locus + " " + ref + " " + alt + " " + str(pval) + " " + str(OR_rr_ra) + " " + str(OR_rr_aa) + " " + str(OR_ra_aa) + "\n")
    print(" OK!")
    if args.manhattan:
        multitools.manhattan(association_filename)
    if args.qqplot:
        multitools.qqplot(association_filename)

if not args.assoc and (args.manhattan or args.qqplot):
    print('Cannot draw plot(s) without performing association test.')
        
if args.info:
    print("Initiating get_info...")
    snp_code = args.info
    print(getinfo.get_info(snp_code, controls_to_analyze))  
               
end = timeit.default_timer()
runtime = end - start
print(runtime)
