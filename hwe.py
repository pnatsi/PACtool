import multitools
from scipy import stats

def hwe(controls, cases):
    snps_hwe = []
    #FIRSTLY UNIFY CONTROL+CASES DATASETS
    for cl in controls:
        for cs in cases:
            if cl[0] == cs[0]:
                for entry in cs[5:]:
                    cl.append(entry)
                snps_hwe.append(cl)
    #THEN CALCULATE HWE STATISTIC FOR EACH SNP
    pairs = []
    result = []
    for snp in snps_hwe:
        pairs.append([snp[0], multitools.get_string(snp)])      #create pairs [snp_code, string]
    for pair in pairs:
        
        string = pair[1]
        obs_rr = string.count('(')      #observed homozygous ref
        obs_ra = string.count('.')      #observed heterozygous 
        obs_aa = string.count(')')      #observed homozygous alt 
        
        p = (2*obs_rr + obs_ra)/(2*(obs_rr + obs_ra + obs_aa))  #frequency of ref
        q = 1-p                                                 #frequency of alt
        n = len(string)                                         #total positions
        
        exp_rr = p**2 * n              #expected homozygous ref
        exp_ra = 2 * p * q * n         #expected heterozygous
        exp_aa = q**2 * n              #expected homozygous alt
        
        s = stats.chisquare([obs_rr, obs_ra, obs_aa], [exp_rr, exp_ra, exp_aa])
        result.append([pair[0], s[0], s[1]])
    return result
    
