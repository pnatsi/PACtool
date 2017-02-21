#    For a given individual that has "B" compute the odds that the same individual has "A"
#    For a given individual that does not have "B" compute the odds that the same individual has "A"
#    Divide the odds from step 1 by the odds from step 2 to obtain the odds ratio (OR).

#    In a more technical language, the OR is a measure of effect size, 
#    describing the strength of association or non-independence between
#    two binary data values

import numpy as np
import multitools
import scipy.stats


def genotype_test(snpCtrl,snpCss): #
    control = multitools.get_string(snpCtrl)
    case = multitools.get_string(snpCss)
    a1 = snpCtrl[3] #allele 1
    a2 = snpCtrl[4] #allele 2
    
    gM = np.array([[case.count('('), case.count('.'), case.count(')')],
                    [control.count('('),control.count('.'),control.count(')')]])
    #gM (genotype Matrix) includes all 3 possible genotypes in two rows
    #one for case and one for control
    # ( = 100, . = 010, ) = 001
    
    OddsRatio={}
    try:
        OddsRatio[a1+a1+'/'+a1+a2]=(gM[0][0]*gM[1][1])/(gM[1][0]*gM[0][1]) #for example TT/TC
        OddsRatio[a1+a1+'/'+a2+a2]=(gM[0][0]*gM[1][2])/(gM[1][0]*gM[0][2])
        OddsRatio[a1+a2+'/'+a2+a2]=(gM[0][1]*gM[1][2])/(gM[1][1]*gM[0][2])
    except RuntimeWarning:
        print('This snp could not be included in association test: ' + snpCtrl[0])    
    try:
        pval = scipy.stats.chi2_contingency(gM)[1]
    except ValueError:
        pval = 'NaN' #there are cases were a genotype will not appear at all so the calculation of p value is impossible
    return(pval, OddsRatio)
