import sys
import math
import concat_ref
import re
import csv
import demux
import mappy as mp
import datetime
import argparse

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", "-i", dest="ref", help="sample fasta reference files", nargs='+', metavar="FILE")
    parser.add_argument("-q", "--query", dest="query", help="fastq long reads file to be sorted", metavar="FILE")
    parser.add_argument("-x", dest="preset", help="preset: sr, map-pb, map-ont, asm5, asm10 or splice")
    parser.add_argument("-n", type=int, dest="min_cnt", help="mininum number of minimizers")
    parser.add_argument("-m", type=int, dest="min_sc", help="minimum chaining score")
    parser.add_argument("-g", type=int, dest="g", help="average genome size")
    parser.add_argument("-r", type=int, dest="bw", help="band width")
    parser.add_argument("-t", type=int, dest="maq", default = 1, help="minimum quality score")
    parser.add_argument("-u", type=int, dest="u", default = 20, help="quality difference to qualify as unique")
    parser.add_argument("-c", dest = "consec", help="concat_ref, ids, fasta# to skip reference concatenation for consecutive runs")
    parser.add_argument("remainder", nargs=argparse.REMAINDER)


    options = parser.parse_args()

    if (not options.ref) or (not options.query):   # if filename is not given
        parser.error('reference.fasta and query.fastq parameters not provided (-i and -q)')
        sys.exit(1)
        
    references = []
    for i in options.ref:
        references.append(i)
    
    if options.remainder:
        print('Unused arguments :' + options.remainder)
    
    print("\nStart: " + str(datetime.datetime.now()))    
    

    k,w = None, None
    if options.g:
        k = math.ceil(math.log((options.g*100),4))
        w = k
        
    if not options.consec:
    #concatenate reference files to single file
        cat_ref, ids, fasta_total = concat_ref.main(references)
    else:
        catlist = options.consec.split(',')
        cat_ref = catlist[0]
        ids = catlist[1]
        fasta_total = catlist[2]
        
    print("\nAligning: " + str(datetime.datetime.now()))
        
    a = mp.Aligner(cat_ref, preset=options.preset, min_cnt=options.min_cnt, min_chain_score=options.min_sc, k=k, w=w, bw=options.bw)
    if not a: raise Exception("ERROR: failed to load/build index file '{}'".format(options.ref))
        
    with open('mapped_seqs.txt', 'w') as output:
	list_of_lists = []
	for name, seq, qual in mp.fastx_read(options.query):    # read one sequence
	    for h in a.map(seq): # traverse hits
            #aligns seq against the index
            	list_of_lists.append(re.split(r'\t+',(name+"\t"+str(h)).rstrip('\t')))
	wr = csv.writer(output)
    	wr.writerows(list_of_lists)
    demux.main('mapped_seqs.txt', fastq_file = options.query, ids=ids, read_count = fasta_total, maq = options.maq, uniqueness = options.u)
    print("\nEnd: " + str(datetime.datetime.now()))
    
    if not options.consec:
        return cat_ref, ids, fasta_total
            
if __name__ == "__main__":
	main(sys.argv)
