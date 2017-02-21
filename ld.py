'''
LD is a measurement of how associated two loci are, i.e. when their frequencies
in a population are somehow dependent from each other.

There are different ways to express LD, two of them being D' and r-squared.

We use the EM algorithm to estimate the frequency of the 'ref-ref' haplotype
and then we use this estimate to calculate D, D' and r-squared.

The initial arbitrary value for p was set to p_r1 * p_r2 (reference haplotypes in two query SNPs), which is the biologically
expected value at someone with equally random haplotype.
'''


import multitools
    
def ld(snp1, snp2, controls, cases):
    #FIND QUERY SNPS IN CONTROLS AND CASES DATASETS
    snp1_cl = [x for x in controls if x[0] == snp1]
    snp2_cl = [x for x in controls if x[0] == snp2]
    snp1_cs = [x for x in cases if x[0] == snp1]
    snp2_cs = [x for x in cases if x[0] == snp2]
    #UNIFY THEM INTO ONE LIST PER SNP
    snp1_total = snp1_cl + snp1_cs[5:]
    snp2_total = snp2_cl + snp2_cs[5:] 
    snp1_final = snp1_total[0] #we dont want a list within a list we 
    snp2_final = snp2_total[0] #just want the one list.
    #CREATE CHARACTERISTIC STRING TO CALCULATE GENOTYPE COUNTS
    pair_snp1 = [snp1_final[0], multitools.get_string(snp1_final)]
    pair_snp2 = [snp2_final[0], multitools.get_string(snp2_final)]
    #INITIALIZE GENOTYPE COUNTS
    n_rrrr = 0
    n_rrra = 0
    n_rraa = 0
    n_rarr = 0
    n_rara = 0
    n_raaa = 0
    n_aarr = 0
    n_aara = 0
    n_aaaa = 0
    #CALCULATE GENOTYPE COUNTS
    string1 = pair_snp1[1]
    string2 = pair_snp2[1]
    for i in range(0, len(string1)):
        if string1[i] == '(':
            if string2[i] == '(':
                n_rrrr += 1
                continue
            elif string2[i] == '.':
                n_rrra += 1
                continue
            elif string2[i] == ')':
                n_rraa += 1
                continue
        elif string1[i] == '.':
            if string2[i] == '(':
                n_rarr += 1
                continue
            elif string2[i] == '.':
                n_rara += 1
                continue
            elif string2[i] == ')':
                n_raaa += 1
                continue
        elif string1[i] == ')':
            if string2[i] == '(':
                n_aarr += 1
                continue
            elif string2[i] == '.':
                n_aara += 1
                continue
            elif string2[i] == ')':
                n_aaaa += 1
                continue
    #TOTAL        
    n = n_rrrr + n_rrra + n_rraa + n_rarr + n_rara + n_raaa + n_aarr + n_aara + n_aaaa
    #CALCULATE SINGLE LOCUS HAPLOTYPES
    p_r1 = (n_rrrr + n_rrra + n_rraa + (n_rarr + n_rara + n_raaa)/2) / n
    p_a1 = 1 - p_r1
    p_r2 = (n_rrrr + n_rarr + n_aarr + (n_rrra + n_rara + n_aara)/2) / n
    p_a2 = 1 - p_r2 
    #SET ARBITRARY P_RR
    p_rr = p_r1 * p_r2
    #EM ALGORITHM
    counter = 0
    while True:
        counter += 1
        if counter > 10000:
            print("EM algorithm failed to converge after 10.000 loops, receiving current output...")
            break
        p_old = p_rr
        try:
            E = 2*n_rrrr + n_rrra + n_rarr + (p_rr*(1 + p_rr - p_r1 - p_r2) * n_rara)/((p_r1 - p_rr)*(p_r2 - p_rr) + p_rr*(1 + p_rr - p_r1 - p_r2))
            p_new = E/(2*n)
        except ZeroDivisionError:
            print("there was an error in the input")
            break
        if abs(p_new - p_old) < 0.0000001:
            break
        else:
            p_rr = p_new
    #CALCULATE D
    D = p_rr - p_r1*p_r2
    #CALCULATE D'
    if D >= 0:
        d_max = min([p_r1*p_a2, p_a1*p_r2])
    else:
        d_max = min([p_r1*p_r2, (1 - p_r1)*(1 - p_r2)])
    try:
        D_prime = D / d_max  
    except ZeroDivisionError:
        D_prime = "D' could not be calculated"
    #CALCULATE R-SQUARED   
    try:
        r_sq = (D**2) / p_r1*p_a1*p_r2*p_a2
    except ZeroDivisionError:
        r_sq = "R-squared could not be calculated"
        
    return [D_prime, r_sq]
