# simulate reads
#!/usr/bin/Rscript

# args <- commandArgs(trailingOnly = TRUE)
# 
# n_args <- length(args)
# 
# read_len <- args[1]
# n_reads <- args[2]
# files <- args[3:n_args]

library(seqinr)
library(ShortRead)
library(Biostrings)

insert <- function(string, idx, ins){
  front <- substr(string, start = 1, stop = idx)
  back <- substr(string, start = (idx+1), stop = nchar(string))
  new_string <- paste(front, ins, back, sep = '')
  return(new_string)
}

delete <- function(string, idx){
  front <- substr(string, start = 1, stop = (idx-1))
  back <- substr(string, start = (idx+1), stop = nchar(string))
  new_string <- paste(front, back, sep = '')
  return(new_string)
}

substitute <- function(string, idx, sub){
  front <- substr(string, start = 1, stop = (idx-1))
  back <- substr(string, start = (idx+1), stop = nchar(string))
  new_string <- paste(front, sub, back, sep = '')
  return(new_string)
}

sample_read <- function(genome, read_len){
  genome_len <- nchar(genome)
  read_len <- read_len - 1
  start <- sample(1:(genome_len - read_len),1)
  read <- substr(genome, start, start+read_len)
  return(read)
}

read_error <- function(read, ins_prob = 0.0002, sub_prob = 0.0002, del_prob = 0.001, prob_sample_size = 10^6){
  nucleotides <- c('a', 't', 'c', 'g')
  
  read_len <- nchar(read)
  error_sample <- sample(1:prob_sample_size, read_len)
  ins_idxs <- which(error_sample < (ins_prob*prob_sample_size))
  error_sample <- sample(1:prob_sample_size, read_len)
  sub_idxs <- which(error_sample < (sub_prob*prob_sample_size))
  error_sample <- sample(1:prob_sample_size, read_len)
  del_idxs <- which(error_sample < (del_prob*prob_sample_size))
  
  for (i in ins_idxs){
    ins <- sample(nucleotides, 1)
    read <- insert(read, i, ins)
  }
  for (i in sub_idxs){
    sub <- sample(nucleotides, 1)
    read <- substitute(read, i, sub)
  }
  for (i in del_idxs){
    read <- delete(read, i)
  }
  
  return(read)
}  

draw_genome_reads <- function(genome, mean_read_length = 1500, num_reads = 10^6){
  genome_len <- nchar(genome)
  read_length_samples <- rpois(num_reads, mean_read_length)
  read_length_samples[which(read_length_samples > genome_len)] <- genome_len
  
  reads <- array(data = NA, dim = c(num_reads))
  q_scores <- array(data = NA, dim = c(num_reads))
  
  for (i in 1:num_reads){
    og_read <- sample_read(genome, read_len = read_length_samples[i])
    # print(grepl(og_read, genome, fixed = TRUE))
    read <- read_error(og_read, ins_prob = 0, sub_prob = 0, del_prob = 0)
    q <- paste(rep('Z', nchar(read)), collapse = '')
    # if (!identical(og_read, read)){
    #   print('error')
    #   print(paste(og_read, read))
    # }
    
    reads[i] <- og_read
    q_scores[i] <- q
  }
  
  list(reads = reads, length = read_length_samples, q_scores = q_scores)
}

sample_contigs_from_genome <- function(genome, mean_contig_length  = 100000, num_contigs = 1, genome_name = 'test'){
  contigs <- draw_genome_reads(genome, mean_read_length = mean_contig_length, num_reads = num_contigs)$reads
  write.fasta(contigs, file.out = paste('sample_contig_', genome_name, collapse = ''), names = 'NA')
  return(contigs)
}

fasta_to_string <- function(file){
  genome <- read.fasta(file)
  # characters <- genome[[1]]
  characters <- c()
  for(i in 1:length(genome)){
    characters <- c(characters, genome[[i]])
  }
  return(paste(characters, collapse = ''))
}

fastq_to_string <- function(file){
  genome <- readFastq(file)@sread
  # characters <- genome[[1]]
  characters <- c()
  for(i in 1:length(genome)){
    read <- toString(genome[[i]])
    characters <- c(characters, read)
  }
  return(paste(characters, collapse = ''))
}

write_fastq <- function(reads, qualitySet, output_file = 'output.fastq', name = 'na'){
  if (length(reads) != length(qualitySet)){return()}
  for (i in 1:length(reads)){
    read <- toupper(reads[i])
    # print(read)
    quality <- qualitySet[i]
    name_line <- paste('@', name, sep = '')
    # quality_line <- paste(rep(quality_char, read_len), sep = '')
    
    write(x = name_line, file = output_file, append = TRUE)
    write(x = read, file = output_file, append = TRUE)
    write(x = '+', file = output_file, append = TRUE)
    write(x = quality, file = output_file, append = TRUE)
  }
}

multiplex_genomes <- function(files, mean_read_length = 2500, num_reads = 10^6, matching_file = 'file_matching.csv'){
  
  mixed_reads <- c()
  mixed_qscores <- c()
  
  for (file in files){
    genome <- fasta_to_string(file)
    write.fasta(genome, file.out = paste('source_', file, sep = ''), names = 'NA')
    # contig <- substr(genome, 1, 100000) #sample_contigs_from_genome(genome = genome, genome_name = file)
    # write.fasta(contig, file.out = paste('sample_contig_', file, collapse = ''), names = 'NA')
    sample_reads <- draw_genome_reads(genome, mean_read_length = mean_read_length, num_reads = num_reads)
    mixed_reads <- c(mixed_reads, sample_reads$reads)
    mixed_qscores <- c(mixed_qscores, sample_reads$q_scores)
    file_name <- paste('sampled', file, sep = '_')
    fastq_file <- strsplit(file, split = 'fasta')[[1]][1]
    fastq_file <- paste('source_', fastq_file, 'fastq', sep = '')
    line <- paste(file_name, fastq_file, sep = ', ')
    write(x = line, file = matching_file, append = TRUE)
    file_matching <- c(matching_file, list(file_name, fastq_file))
    write.fasta(as.list(sample_reads$reads), file.out = file_name, names = 'NA', open = 'w')
  }
  
  # for (i in 1:length(file_matching)){
  #   files <- file_matching[1]
  #   
  #   
  # }
  
  shuffle_reads <- sample(1:length(mixed_reads), size = length(mixed_reads))
  
  mixed_reads <- mixed_reads[shuffle_reads] #sample(mixed_reads, size = length(mixed_reads))
  mixed_qscores <- mixed_qscores[shuffle_reads]
  ids <- rep('name', length(mixed_reads))
  
  write_fastq(reads = mixed_reads, qualitySet = mixed_qscores)
  
  # 
  # dna <- DNAStringSet(mixed_reads)
  # qual <- BStringSet(mixed_qscores)
  # ids <- BStringSet(ids)
  # 
  # reads <- ShortReadQ(sread = dna, quality = qual, id = ids)
  # writeFastq(reads, file = 'multiplexed_reads.fastq')
  # write.fasta(as.list(mixed_reads), file.out = 'multiplexed_reads.fasta', names = 'NA')
}

# genome <- fasta_to_string('sobrinus_NCTC10921.fasta')
# files <- c('pneumoniae_R6.fasta', 'pyogenes_M1_GAS.fasta', 'sobrinus_NCTC10921.fasta')
files <- c('10919contigs.fasta', '6715_15contigs.fasta', 'SL1contigs.fasta')
# genome <- fasta_to_string(files[1])
# sample_reads <- draw_genome_reads(genome, mean_read_length = 10, num_reads = 100)
multiplex_genomes(files, mean_read_length = 6000, num_reads = 10)
