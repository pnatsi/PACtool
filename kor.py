'''
We define four functions for pre-analysis input file manipulation.
With these keep/remove snps and keep/remove samples actions
can be applied to the input file.
'''

def keep_snps(gen_file, keep_file):
    gen = []
    keep = []
    pop=[]
    with open(gen_file, 'r') as g:
        for line in g.readlines():
            gen.append(line.strip())
    with open(keep_file, 'r') as k:
        for line in k.readlines():
            keep.append(line.strip())
    for snp in gen:
        splitted = snp.split()
        if splitted[0] in keep:             #we check if the snp_code of the initial file
            pop.append(snp)                 #is in the keepsnps file and we obtain it
    return pop
    
def remove_snps(gen_file, remove_file):
    gen = []
    remove = []
    with open(gen_file, 'r') as g:
        for line in g.readlines():
            gen.append(line.strip())
    with open(remove_file, 'r') as r:
        for line in r.readlines():
            remove.append(line.strip())
    for snp in gen:
        splitted = snp.split()              #we check if the snp_code of the initial file
        if splitted[0] in remove:           #is in the removesnps file and we remove it
            gen.remove(snp)
    return gen

def keep_samples(clgen_file, csgen_file, keeps_file):
    clgen = []
    csgen = []
    keeps = []
    clpop = []
    cspop = []
    i = 5
    #first we read the input files and append their lines to corresponding lists
    with open(clgen_file, 'r') as g:            
        for line in g.readlines():
            clgen.append(line.strip())
    with open(csgen_file, 'r') as g:
        for line in g.readlines():
            csgen.append(line.strip())        
    with open(keeps_file, 'r') as k:
        for line in k.readlines():
            keeps.append(line.strip())
    #then we isolate the control/case codes that want to be kept
    cl = [x[8:] for x in keeps if x[0:7] == 'control']
    cs = [x[5:] for x in keeps if x[0:4] == 'case']
    #and here we obtain the final snps with only the preferred samples
    for snp in clgen:
        splitted = snp.split()                                  #split the line jnto list
        topop = splitted[:i]                                    #first 5 elements go as they were
        for number in cl:                                       #and with this for loop we keep
            topop.append(splitted[i+(int(number)-1)*3])         #only the preferred samples(triplets)    
            topop.append(splitted[i+(int(number)-1)*3+1])
            topop.append(splitted[i+(int(number)-1)*3+2])
        clpop.append(topop)
    for snp in csgen:                                           #same thing for cases
        splitted = snp.split()
        topop = splitted[:i]
        for number in cs:
            topop.append(splitted[i+(int(number)-1)*3])
            topop.append(splitted[i+(int(number)-1)*3+1])
            topop.append(splitted[i+(int(number)-1)*3+2])
        cspop.append(topop)
    return clpop, cspop

def remove_samples(clgen_file, csgen_file, removes_file):
    clgen = []
    csgen = []
    removes = []
    clpop = []
    cspop = []
    i = 5
    #first we read the input files and append their lines to corresponding lists
    with open(clgen_file, 'r') as g:
        for line in g.readlines():
            clgen.append(line.strip())
    with open(csgen_file, 'r') as g:
        for line in g.readlines():
            csgen.append(line.strip())        
    with open(removes_file, 'r') as r:
        for line in r.readlines():
            removes.append(line.strip()) 
    #then we isolate the control/case codes that want to be kept
    cl = [x[8:] for x in removes if x[0:7] == 'control']
    cs = [x[5:] for x in removes if x[0:4] == 'case']
    #and here we obtain the final snps by removing the unpreferred samples
    for snp in clgen:
        splitted = snp.split()                                          #splitting the line into list
        for number in cl:                                               
            del splitted[i+(int(number)-1)*3:i+(int(number)-1)*3+3]     #here we delete the triplets that
            clpop.append(splitted)                                      #indicated in the removesamples file
    #same thing for cases        
    for snp in csgen:
        splitted = snp.split()
        for number in cs:
            del splitted[i+(int(number)-1)*3:i+(int(number)-1)*3+3]
            cspop.append(splitted) 
    #a final action to remove possible duplicates from the output
    clfinal = []
    for i in clpop:
        if i not in clfinal:
            clfinal.append(i)
    csfinal = []
    for i in cspop:
        if i not in csfinal:
            csfinal.append(i)
    return clfinal, csfinal

    