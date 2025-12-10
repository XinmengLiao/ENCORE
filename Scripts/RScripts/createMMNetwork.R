library(dplyr)
library(tidyr)

sys.args <- commandArgs(trailingOnly = TRUE)
taxo_dir <- sys.args[1]
fva_dir <- sys.args[2]
output_file <- sys.args[3]

# list taxonomy files
taxo.all <- data.frame()
taxo_files <- list.files(taxo_dir)
print("Now reading taxonomy files from all samples. ")
for (i in taxo_files){
  sampleid = unlist(strsplit(i,"_gtdbtk_summary"))[1]
  taxo.tmp <- read.csv(file.path(taxo_dir,i),header = T,sep = "\t") %>% 
    mutate(Sample = sampleid) %>% unique()
  taxo.all <- rbind(taxo.all, taxo.tmp) 
}
taxo.all <- taxo.all %>% select(Sample, everything()) %>% 
  mutate(user_genome = gsub("\\.","",user_genome)) %>% 
  rename(Bin = user_genome) %>% 
  filter(!is.na(classification)) 

# list fva files
fva.all <- data.frame()
fva_files <- list.files(fva_dir)
print("Now reading FVA files from all samples. ")
for (i in fva_files){
  sampleid = unlist(strsplit(i,"\\."))[1]
  binid = gsub("\\FVA","", unlist(strsplit(i,"\\."))[2])
  fva.tmp <- read.csv(file.path(fva_dir,i),header = T,sep = "\t") %>% 
    mutate(Sample = sampleid, Bin = binid) %>% unique() 
  fva.all <- rbind(fva.all, fva.tmp) 
}

# filter only mandatory exchange reactions 
fva.all.exchange <- fva.all %>% 
  filter(UB_exc < 0 | LB_exc > 0 ) %>% unique() %>% 
  mutate(Direction = case_when(UB_exc < 0 ~ "Uptake",
                               LB_exc > 0 ~ "Secretion",
                               TRUE ~ "None"))

write.table(fva.all, file = file.path(output_file, "fva_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fva.all.exchange, file = file.path(output_file, "fva_all_nonzero.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# merge taxonomy and fva table together

merged_taxa_fva <- fva.all.exchange %>% 
  left_join(., taxo.all %>% select(Sample, Bin,classification) %>% unique(), by  = c("Sample","Bin")) %>% 
  separate(classification,
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";",
           remove = FALSE) %>%
  mutate(across(Domain:Species, ~gsub("^[a-z]__", "", .))) %>% 
  filter(!is.na(Species)) %>% unique()

create.network <- function(data){
  secrete <- data %>% filter(Direction == "Secretion") %>%
    select(classification,exchangeRxns,Direction)  %>% 
    unique() %>% 
    dplyr::rename(Vertex.A = classification, Vertex.B = exchangeRxns) %>% 
    mutate(Type.A = "Microbe", Type.B = "Metabolite")
  network <- data %>% filter(Direction == "Uptake") %>% 
    select(exchangeRxns,classification,Direction) %>% 
    unique() %>% 
    dplyr::rename(Vertex.A = exchangeRxns, Vertex.B = classification) %>% 
    mutate(Type.A = "Metabolite", Type.B = "Microbe") %>% 
    rbind(., secrete) %>% 
    unique() %>% 
    filter(Vertex.A != Vertex.B) %>% 
    group_by(edge_group = paste(pmin(Vertex.A, Vertex.B), pmax(Vertex.A, Vertex.B))) %>% 
    filter(n() == 1) %>%  
    ungroup() %>% 
    select(Vertex.A, Type.A, Vertex.B, Type.B,everything()) %>% select(-edge_group)
  return(network)
}

network.all <- create.network(merged_taxa_fva)
write.table(network.all, file = file.path(output_file, "MMNetwork.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

