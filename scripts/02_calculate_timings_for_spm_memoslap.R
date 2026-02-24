
############################################################
# Script: Extract Block Timings for SPM
#
# Author: Gerasimos Gerardos
# Project: MeMoSLAP Project (P04–Leipzig)
# Paper: Task load modulates network interactions between
#        bilateral fronto-parietal and cerebellar areas
#        during verbal working memory
#
# Description:
# This script extracts block onset times and block durations
# from behavioral task CSV files and formats them for SPM
# first-level fMRI analyses. It processes multiple participant
# files, aligns timing information to the experiment start,
# computes onset times and durations for task, fixation, and
# instruction blocks, and exports both participant-level and
# group-level timing tables.
#
# Output:
# - Individual block onset files (*_block_starts.csv)
# - Individual duration files (*_block_durations.csv)
# - Combined files (all_starts.csv, all_durations.csv)
#
# Date: 2026-02-11
############################################################

work_directory <- "/data/behavior/long_task_table/" # change to your directory 

# if you chose a directory, the working directory will be set to this path and the text informs that it worked
if (nzchar(work_directory)) {
  cat("Selected directory:", work_directory, "\n")
  setwd(work_directory)
}

L = list.files(work_directory, pattern = "*.csv")

for (j in 1:length(L)) {
  
  data_csv <- read.csv2(L[j], head = TRUE, sep = ",", na = "NA", stringsAsFactors = FALSE)
  
  
  start_exp <-  as.numeric(data_csv[5, 37]) + as.numeric(data_csv[5, 36])

  timings <- data_csv[6:nrow(data_csv), c('block.thisRepN', 'trials.thisRepN', 'n_trial_tot', 'curr_n_back', 'letter.started', 
                                              'letter.stopped', 'void.started', 'void.stopped', 'cross.started', 'text.started', 'ibi.duration')]
  for (i in 1:ncol(timings)) {
    timings[ , i] <- as.numeric(timings[ , i],  na.rm = TRUE)
  }
  

  block_starts <- matrix(nrow = 12, ncol = 8)
  # rownames(block_starts) <- c("0back", "1back", "2back", "3back", "0fix", "1fix", "2fix", "3fix", "0instr", "1instr", "2instr", "3instr")
  colnames(block_starts) <- c("1", "2", "3", "4", "5", "6", "block", "session")
  block_starts[,c("session")] <- data_csv[1,c("participant")]
  block_starts[,c("block")] <- c("0back", "1back", "2back", "3back", "0fix", "1fix", "2fix", "3fix", "0instr", "1instr", "2instr", "3instr")
  
  durations <- matrix(nrow = 12, ncol = 8)
  # rownames(durations) <- c("0back", "1back", "2back", "3back", "0fix", "1fix", "2fix", "3fix", "0instr", "1instr", "2instr", "3instr")
  colnames(durations) <- c("1", "2", "3", "4", "5", "6", "block", "session")
  durations[,c("session")] <- data_csv[1,c("participant")]
  durations[,c("block")] <- c("0back", "1back", "2back", "3back", "0fix", "1fix", "2fix", "3fix", "0instr", "1instr", "2instr", "3instr")
  
  #outcomes <- matrix(nrow = 16, ncol = length(timings[,1]))
  
  for (i in 0 : 3) {
    start_t = as.numeric(timings[timings$curr_n_back==i & timings$trials.thisRepN==0 & !is.na(timings$trials.thisRepN), c("letter.started")]) - start_exp
    block_starts[i+1, 1:length(start_t)] = start_t
    
    start_in = as.numeric(timings[timings$curr_n_back==i & timings$trials.thisRepN==0 & !is.na(timings$trials.thisRepN), c("text.started")]) - start_exp
    block_starts[i+9, 1:length(start_in)] = start_in
    
    start_fix = as.numeric(timings[which(timings$curr_n_back==i & timings$trials.thisRepN==34)+1, c("cross.started")]) - start_exp
    block_starts[i+5, 1:length(start_in)] = start_fix
    
    end_t <- as.numeric(timings[timings$curr_n_back==i & timings$trials.thisRepN==34 & !is.na(timings$trials.thisRepN), c("void.stopped")]) - start_exp 
    
    nduration <- end_t - start_t 
    nduration[is.na(nduration)] <- mean(nduration[!is.na(nduration)])
    durations[i+1, 1:length(nduration)] = nduration
    
    iduration = start_t - start_in
    durations[i+9, 1:length(iduration)] = iduration
    
    start_plus1 <- as.numeric(timings[which(timings$curr_n_back==i & timings$trials.thisRepN==34)+1, c("text.started")]) - start_exp
    
    fduration <- start_plus1
    fduration[!is.na(start_plus1)] = start_plus1[!is.na(start_plus1)] - start_fix[!is.na(start_plus1)]
    fduration[is.na(start_plus1)] = 1420 - start_fix[is.na(start_plus1)]
    durations[i+5, 1:length(fduration)] = fduration
  
  }
  
  start_path <- paste("/results/fMRI_analysis/timings_SPM/block_starts/", data_csv[1,c("participant")], "_block_starts.csv", sep = "") #change to your directory 
  dur_path <- paste("/results/fMRI_analysis/timings_SPM/block_durations/", data_csv[1,c("participant")], "_block_durations.csv", sep = "") #change to your directory 
  write.csv(block_starts, file = start_path, row.names = FALSE)
  write.csv(durations, file = dur_path, row.names = FALSE)
  
  
  if (j == 1) {
    all_starts <- block_starts
    all_durations <- durations
  } else {
    all_starts <- rbind(all_starts, block_starts)
    all_durations <- rbind(all_durations, durations)
  }
}

start_path <- "/results/fMRI_analysis/timings_SPM/block_starts/all_starts.csv" #change to your directory 
dur_path <- "/results/fMRI_analysis/timings_SPM/block_durations/all_durations.csv" #change to your directory 
write.csv(all_starts, file = start_path, row.names = FALSE)
write.csv(all_durations, file = dur_path, row.names = FALSE)

