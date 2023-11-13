library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(janitor)

mydir <- getwd()
print(mydir)


# Upload file
args <- commandArgs(trailingOnly = TRUE)
excel_file <- args[1]
my_data <- read_excel(excel_file)

# Filter: complete genomes
filtered_data <- function(data_name){
  data = data_name %>% 
    clean_names() %>%
    mutate(genome_id = as.character(genome_id)) %>%
    filter((genome_status == 'Complete') %>% replace_na(TRUE)) %>% 
    filter((genome_quality == 'Good') %>% replace_na(TRUE)) %>%
    mutate(plasmids = replace_na(plasmids,0)) %>% 
    mutate(nm_contigs = contigs - plasmids) %>% 
    filter(nm_contigs <= 100)
  nm_sp <- nrow(data)
  print(nm_sp)
  return(data.frame(data))
}

#filtered <- lapply(my_data, filtered_data)
filtered_complete <- filtered_data(my_data)


#Write csv and tsv file
write.csv(filtered_complete, "filtered_complete.csv")
print("File created")

data_c <- filtered_complete %>% 
  select(genome_id) %>% 
  mutate(genome_id = paste0(genome_id, ".fna"))

write_delim(data_c, "genome_list.txt", delim = "\t", col_names = FALSE)

#Filter WGS genomes
wgs_data <- function(data_name,i){
  data = data_name %>% 
    clean_names() %>% 
    filter((check_m_completeness > 90) %>% replace_na(TRUE)) %>% 
    filter((check_m_contamination < 10) %>% replace_na(TRUE)) %>%
    mutate(plasmids = replace_na(plasmids,0)) %>% 
    mutate(nm_contigs = contigs - plasmids) %>% 
    filter(nm_contigs <= 100) 
  nm_sp <- nrow(data)
  print(nm_sp)
  return(data.frame(data))
}

#wgs <- lapply(my_data, wgs_data)
filtered_wgs <- wgs_data(my_data)

#Write csv
write.csv(filtered_wgs, "filtered_wgs.csv")

#Write tsv file - use it for download .fna files from BV-BRC database
data_w <- filtered_wgs %>% 
  select(genome_id) %>%
  mutate(genome_id = paste0(genome_id, ".fna"))

write_delim(data_w, "Data/genome_list_wgs.txt", delim = "\t", col_names = FALSE)