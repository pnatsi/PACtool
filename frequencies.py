import timeit

start = timeit.default_timer()

def frequencies(snp):

    code = snp[0]

    string = multitools.get_string(snp)    
    ref_freq = (string.count('(')/len(string) + (string.count('.')/len(string)))/2
    alt_freq = (string.count(')')/len(string) + (string.count('.')/len(string)))/2
    final = [code, ref_freq, alt_freq]
    return final
