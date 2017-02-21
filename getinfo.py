import requests
from pyliftover import LiftOver
import sys

def get_info(snp_code, snp_list):
    for snp in snp_list:
        if snp[0] == snp_code:
            hg18_coordinates = snp[2]
            ref = snp[3]
            alt = snp[4]
            pylift_id = int(hg18_coordinates) - 1
            lo = LiftOver('hg18', 'hg38')             
            sys.stdout.write("\n\nUpdating chromosome position from ncbi build 36(hg18) to build 38(hg38)..")                
            pylift_tuple = lo.convert_coordinate('chr20', pylift_id)
                 
            sys.stdout.write("  ✔\n")
                 
            hg38_coordinates = int(pylift_tuple[0][1]) + 1
                                        
            chr_prefix = '20:g.'
            request_id = chr_prefix+str(hg38_coordinates)+ref+'>'+alt+'?'
            url = 'http://rest.ensembl.org/vep/human/hgvs/'
            headers={ "Content-Type" : "application/json"}
                 
            r = requests.get(url+request_id, headers=headers)
            if (r.ok) == True:
                data = r.json()
                print ('\nQuery successful!\n\nYou can find info about', snp_code, 'below:\n_____________________________________\n')
                print ('\nYour json contains the following', len(data[0].keys()) , 'keys:\n')
                print ('\n'.join(data[0].keys()), '\n\n')
                print (data)
            else:
                print ('Sorry, information for ', snp_code, " currently unavailable.")
            break
        else:
            print('Sorry, but', snp_code, 'is not included in the dataset.\nOr maybe you have misspelled the snp_id?.\nPlease try again.')
    return ""

        
        