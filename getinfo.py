#!/usr/bin/env python3

# THE FOLLOWING PROGRAMME IS PART OF the PACtool package (https://pypi.python.org/pypi/PACtool)
 
# FUNCTIONALITY:
# This programme can be used to automatically fetch information about 
# a selected variant from Ensenmbl's Variant Effect Prredictor database.
# The current version updates by default coordinates from ncbi build36 (hg18) to build38 (hg38)
 
# OUTPUT FILE:
# RETRIEVED INFORMATION IN JSON FORMAT IS SAVED IN A FILE WITH FILENAME: snp_code.info 


import requests
from pyliftover import LiftOver
import sys

# snp_list is by default assigned to a list created in previous steps of PACtool analysis.(controls_to_analyze)


# controls_to_analyze: A LIST with as many sublists as the snps. 
#                      Each list sublist has 5 elements, the info of the first 5 columns of the files gwas.*.gen
#                      For the selected SNPs in the keep/remove steps.


def get_info(snp_code, snp_list):

    
    chr_prefix = "20"              
    old_build  = "hg18"           
    new_build  = "hg38" 
    
    cols=[x for x in zip(*snp_list)] # This bit, unzips [* asterisk] the lis, which is actually a LIST of lists
                                     # and 'dumps' all same indexed elements from sublists in a tuple,
                                     # The list cols, contains 5 tuples; len(cols)=5
                                    
                                    
    
                                     
                                   
    snp_codes        = cols[0]      # snp_codes         :  cols[0] = ('snp_0', 'snp_1', 'snp_2', 'snp_3', ...) __ 1st column of gwas.cases.gen
    hg18_coordinates = cols[2]      # hg18_coordinates  :  cols[2] = ('9098', '9150', '9795', '10731',.......) __ 3rd column of gwas.cases.gen
    ref              = cols[3]      # ref_base          :  cols[3] = ('C', 'T', 'G', 'C', 'A', '.............) __ 4th column of gwas.cases.gen
    alt              = cols[4]      # alt_base          :  cols[4] = ('T', 'A', 'T', 'A', 'C', ..............) __ 5th column of gwas.cases.gen
                                    # the whole rationale behind this, is that the index stays true amongst the wbove lists
                                    # in relation to the initial file.
                                    # So if we took the elements cols[0][0],cols[1][0], cols[2][0], cols[3][0],cols[][0]
                                    # we'd get all the elements in the row of snp_0 in the file gwas.cases.gen                        
                                              

    # NOW THAT WE HAVE EACH COLUMN OF THE FILE STORED IN A TUPLE,
    # WE CAN USE THEM TO ITERATE AND GENERATE THE HGVS id 
    # THAT WE'LL NEED FOR SENDING THE GET REQUEST TO GET INFO.
    
    if snp_codes.count(snp_code)!=0:                                        # this checks if the user_input matches a snp_id i.e. snp_90 from the file 12345_long.txt
        for i in range(0,len(snp_codes)):
            if snp_code == snp_codes[i]:
                
              
                           
                
                sys.stdout.write("\n\nUpdating chromosome coordinates to new assembly..")
                   
                pylift_id = int(hg18_coordinates[i])
                lo = LiftOver(old_build, new_build)         # Stating from ('hg18') to('hg38') which build we want our coordinates to be updated  
                   
                   
                # We'll use: LiftOver.convert_coordinate('chrX', 'XXXXX') to update the position on the chromosome
                # All SNPs are located on chromosome 20, so 
                # the first argument of lo.convert_coordinates will be 'chr20' for all SNPs
                # the second argument of lo.convert_coordinates will be the pylift_id from above
                   
                # The output of lo.convert_coordinates will be a list with 1 element; a tuple with 4 elements
                # i.e.: pylift_tuple = [('chr20', 80456, '+', 5643036713)] 
                # We can access the elements of the tuple by using 2 sq.brackets as index indicators
                # SO:
                # pylift_tuple[0][0] = 'chr20'     ______ chromosome i.e. 
                # pylift_tuple[0][1] = '80456'     ______ coordinates in hg38 <---- we'll need this one
                # pylift_tuple[0][2] = '+'         ______ DNA strand, + for coding, - for non-coding
                # pylift_tuple[0][3] = 5643036713  ______ allignment score*
                   
                # NOTE: We will need the coordinates aka the pylift_tuple[0][1]
                
       
                pylift_tuple = lo.convert_coordinate('chr'+chr_prefix, pylift_id) 
                    
                sys.stdout.write("  Done \n")
                   
            
                # The pylift_tuple[0][1] is a string with the updated coordinates
                # hg38_coordinates, will be used to make the genomic HGVS id, necessary for the request
                   
                hg38_coordinates = int(pylift_tuple[0][1])  #the updated coordinates
                   
                   
                # We have chosen to use GET requests to get info
                # Since NOT all samples have an rs_id 
                # We will try to reconstruct the genomic HGVS id for all samples using their coordinates 
                # Let's take a look at the url link used for a GET request:
                # url = 'http://rest.ensembl.org/vep/human/hgvs/9:g.22125504G>C?'
                # For all of our samples the first part i.e. 'http://rest.ensembl.org/vep/human/hgvs/20:g.
                # would be the same.
                # All we need is for each snp to reconstruct the remaining bit,i.e. 22125504G>C?
                # which is actually a string of this form:
                # i.e. for snp_0: hg38_coordinates + snp[3][0] + '>' + snp[4]+'?'
                   
                             
                   
                # Let request_id 
                # be the variable that holds the reconstructed HGVS id for each SNP                
                request_id = chr_prefix+":g"+str(hg38_coordinates)+ref[i]+'>'+alt[i]+'?'
                url = 'http://rest.ensembl.org/vep/human/hgvs/'
                headers={ "Content-Type" : "application/json"}
                   
                #At last, actually making the request using requests lib:
                r = requests.get(url+request_id, headers=headers) 
                if (r.ok) == True:
                    data = r.json()
                    print ('\nQuery successful!\n\nInfo about', snp_code, "has been saved in the file:  ", snp_code+".info") 
                    print ('\nYour json file contains the following', len(data[0].keys()) , 'keys:\n')
                    print ('\n'.join(data[0].keys()), '\n\n')
     
                    saveout = sys.stdout
                    file = snp_code+".info"
                     
                    save = open(file, 'w')  #saving in a file the output of print which is a decoded json file
                    sys.stdout = save
                     
                    print (data)
                   
                    sys.stdout = saveout
                    save.close()           #closing the file, and now stdout again on terminal
                                      
                else:
                    print ('Sorry, information for ', snp_code, " currently unavailable.")
   
    #REMINDER: The else below, goes with the:
    # if snp[0].count(snp_id)!=0:
    # If user didn't type the right snp_id, 
    # the following message will be printed:
    else:
        print('\nSorry, but', snp_code, 'is not included in the dataset.\nOr maybe you have misspelled the snp_id?.\nPlease try again.\n')
    return("")
     
         
