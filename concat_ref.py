import sys
import numpy
from Bio import SeqIO
from collections import OrderedDict

def main(references):
    lst = []
    mydict = OrderedDict()
    print("\nReferences provided: ")
    for arg in references:
        print(str(arg))
        lst.append(arg)
    for index, item in enumerate(lst):
        mydict['id'+str(index+1)+'_pref'] = str(item) + "_"
  
    ids = []
    
    read_count = 0
    
    with open('all_ref.fasta', 'w') as output:

        for index, label in enumerate(lst):
            fasta_file = SeqIO.parse(label, 'fasta')
            for fasta in fasta_file:
                fasta.id = mydict.values()[index] + fasta.id
                SeqIO.write(fasta, output, "fasta")
                read_count = read_count + 1
            ids.append(mydict.values()[index])
            
    out_file = 'all_ref.fasta'
    return out_file, ids, read_count

if __name__ == "__main__":
    main(sys.argv[1:])  

