# analyze demux
#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

n_args <- length(args)

read_len <- args[1]
n_reads <- args[2]
files <- args[3:n_args]


library(seqinr)
library(ShortRead)
library(Biostrings)

demux_scoring <- function(actual_reads, predicted_reads){
  actual_read_presence <- unlist(lapply(actual_reads, function(x){x %in% predicted_reads@sread}))
  # print(actual_read_presence)
  num_matched <- length(which(actual_read_presence))
  num_missed <- length(which(!actual_read_presence))
  num_mistake <- length(predicted_reads@sread) - num_matched
  
  list(match = num_matched, missed = num_missed, mistake = num_mistake)
}

process_read_files <- function(actual_read_file, predicted_read_file){
  actual_reads <- read.fasta(actual_read_file)
  reads <- c()
  for (i in 1:length(actual_reads)){
    read <- paste(actual_reads[[i]], collapse = '')
    reads <- c(reads, read)
  }
  
  actual_reads <- reads
  predicted_reads <- readFastq(predicted_read_file)
  return(demux_scoring(actual_reads, predicted_reads))
}

calculate_confusion_matrix <- function(files){
  n <- nrow(files)
  
  mtx <- matrix(data = 0, nrow = n, ncol = n)
  rownames(mtx) <- files[,1]
  colnames(mtx) <- files[,2]
  
  for (i in 1:n){
    actual <- files[i,1]
    predicted <- files[i,2]
    
    scores <- process_read_files(actual_read_file = actual, predicted_read_file = predicted)
    
    mtx[i,i] <- scores$match
    
    other <- which(files[,2] != predicted)
    for (j in others){
      scores <- process_read_files(actual_read_file = actual, predicted_read_file = files[j,2])
      mtx[i,j] <- scores$match
    }
  }
  
  return(mtx)
}

matching_file <- 'file_matching.csv'
files <- read.csv(matching_file, header = FALSE)
colnames(files) <- c('actual read file', 'predicted read file')