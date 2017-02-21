import matplotlib.pyplot as plt
import numpy as np 
import scipy.stats as stats

def manhattan(association_output):
    with open(association_output, 'r') as f:
        datlist = [x.split() for x in f]
        datcols = [x for x in zip(*datlist)]
        loci = [int(x) for x in list(datcols[1][1:])]
        pvals = [-np.log10(float(x)) for x in list(datcols[4][1:])]
    
    fig, ax = plt.subplots()
    ax.set_xticks([0,10000000,25000000,40000000,55000000])
    ax.set_xlim(0, 60000000)
    ax.set_ylim(0, 10)
    plt.scatter(loci,pvals)
    plt.ticklabel_format(style = 'plain')
    plt.title("MANHATTAN PLOT")
    plt.xlabel("LOCUS")
    plt.ylabel("-LOG10(P-VALUE)")
    plt.savefig("manhattan.eps")

def qqplot(association_output):
    with open(association_output, 'r') as f:
        datlist = [x.split() for x in f]
        datcols = [x for x in zip(*datlist)]
        pvals = [-np.log10(float(x)) for x in list(datcols[4][1:])]
    
    
    stats.probplot(pvals, dist="norm", plot=plt)
    plt.savefig("qqplot.eps")
    
def get_string(snp):
    i = 5                           #the genotypes start from 6th column
    triplets = []                   #a list with all genotype triplets for the snp

    while i < len(snp):
        triplet = snp[i:i+3]
        triplets.append(str(triplet[0]) + str(triplet[1]) + str(triplet[2]))
        i = i + 3
    
    string = ""                     #this will be the characteristic string for the snp
                                    #to make finding frequencies easier!
    
    for triplet in triplets:        #with this for loop we fill the characteristic string
        if triplet == '100':
            string += '('
        elif triplet == '010':
            string += '.'
        elif triplet == '001':
            string += ')'           #it will look sth like that '((.)(..))...)()'
        else:
            print("There is some error in the columns of the input")
            
    return string
