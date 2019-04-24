import sys
import Bio
from Bio import SeqIO
from operator import itemgetter

def main(mapped, fastq_file, ids, read_count, maq, uniqueness): 
    
    fasta_total = int(read_count)
    l = len(ids)
    files = [0]*l
    for index, label in enumerate(ids):
        newlabel = label.replace('.fasta_','')
        files[index] = newlabel
        
    fastq_count = 0
    #create array of ints to count reads sorted to each file
    label_counts = [0]*l
    #clear all destination files
    for index, label in enumerate(files):
        open(label+'.fastq', 'w').close()
        #open(label+'maq.fastq','w').close()
        open ("unqualified.fastq", 'w').close()

    
    qualified_list = []
    unqualified_list = []
    f = open(mapped, 'r')
    for line in f.readlines():
        parsed_list = line.split(",")
        #mapping quality threshold based on observed mapped_seqs file
        #parsed_list order: name, query_start, query_end, +-strand, ref_node, ref_node_length, ref_start, ref_end, match_length, align_length, mapq, cg:Z:cigar_str
        if int(parsed_list[10]) > int(maq):
            #print(parsed_list)
            qualified_list.append(parsed_list)
        else:
            unqualified_list.append(parsed_list)
            
            
    fastqs = SeqIO.parse(fastq_file, 'fastq')
    approved = []
    
    #write fastqs unused to new file
    unique_unqual = []
    unqual_count = 0
    for item in unqualified_list:
        if (item[0] not in [item2[0] for item2 in qualified_list]) and item[0] not in [item2[0] for item2 in unique_unqual]:
            unique_unqual.append(item)
            unqual_count = unqual_count + 1
    
    for fastq in fastqs:
        #write fastqs unused to new file
        for item in unique_unqual:
            if item[0] == fastq.id:
                with open('unqualified.fastq', 'a') as output:
                    SeqIO.write(fastq, output, "fastq")
        #continue with qualified reads            
        fastq_count = fastq_count+1
        match_count=0
        passed_items = []
        for item in qualified_list:
            if item[0] == fastq.id:
                match_count = match_count+1
                passed_items.append(item)
                
        if match_count > 1:
            #check if one of duplicate matches have high enough relative quality score
            qual_ref = []
            ref = []
            for item in passed_items:
                qual_ref.append([item[4], int(item[10]), item])
                pref = item[4].split(".fasta",1)[0] + ".fasta"
                ref.append(pref)
            qual_ref.sort(key=itemgetter(1), reverse=True)
            if len(list(set(ref)))==1:
                if passed_items[0] not in approved:
                    approved.append(passed_items[0])
                    
            else:
                for index, item in enumerate(qual_ref):
                    if index < len(qual_ref)-1:
                        if item[1] == qual_ref[0][1]:
                            if not (qual_ref[index][0] != qual_ref[index+1][0] and (int(qual_ref[index][1]) < int(qual_ref[index+1][1])+uniqueness)):
                                approved.append(item[2])
            
            for item in approved:
                for index, label in enumerate(files):
                    if label == item[4].split(".fasta")[0]:
                        with open(label+'.fastq', 'a') as output:
                            SeqIO.write(fastq, output, "fastq")
                        label_counts[index] = label_counts[index]+1
                        #with open(label+'maq.fastq', 'a') as output:
                        #    output.write(str(item[0])+ '\t' + str(item[10]) +'\n')
            approved = []
                    
        elif match_count == 1:
                
            #write unique to appropriate file
            for index, label in enumerate(files):
                if label == passed_items[0][4].split(".fasta")[0]:
                    with open(label+'.fastq', 'a') as output:
                        SeqIO.write(fastq, output, "fastq")
                    label_counts[index] = label_counts[index]+1
                    #with open(label+'maq.fastq', 'a') as output:
                    #    output.write(str(passed_items[0][0])+ '\t' + str(passed_items[0][10]) +'\n')
             
    print("\nOutput files containing sorted query reads: ")
    for index, label in enumerate(files):
        print(label + '.fastq' + "\t" + str(label_counts[index]) + "/" + str(fastq_count))
    print("Unqualified reads:\tunqualified.fastq\t" + str(unqual_count) + "/" + str(fastq_count))
        

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
