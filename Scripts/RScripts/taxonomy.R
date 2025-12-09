if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# Taxonomy file from GTDB-tk
files <- list.files("./")
taxo_file <- files[grep(".*summary.tsv", files)]
taxo <- read.csv(taxo_file, header = T,sep = "\t") %>% select(user_genome, classification)
class_colname <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
taxo_split <- as.data.frame(do.call(rbind, strsplit(taxo$classification, split = ";", fixed = TRUE)),stringsAsFactors = F)
colnames(taxo_split) <- class_colname
taxo <- cbind(taxo, taxo_split)
taxo <- taxo %>%  
  mutate(Species = gsub("s__", "", Species),         
			Species = gsub("_E", "", Species),         
			Genus = gsub("g__","",Genus),         
			Genus = gsub("_E","",Genus),
			Species = gsub("_A", "", Species),         
			Genus = gsub("_A","",Genus)) %>% 
  mutate(Species_cap  = toupper(Species)) %>% 
  mutate(Species_cap = gsub("SP.*", "SP.", Species_cap),  
         Species_cap = gsub("_.", "", Species_cap),  )
taxo <- taxo[!duplicated(taxo$user_genome),]

# KEGG database 
kegglist <- read.csv("kegglist.txt",header = T,sep = "\t") %>% 
  mutate(Species_cap = toupper(Speices)) %>% select(1,2,Species_cap)

# Uniprot database
proteome <- read.csv("UniprotID.txt",header = T,sep = "\t") %>% 
  mutate(Brief = gsub("Candidatus","", Brief)) %>% 
  mutate(Brief = gsub(" str.*" , "",Brief)) %>% 
  mutate(Brief = gsub(" subsp.*" , "",Brief)) %>%  
  mutate(Brief = gsub(" serovar.*" , "",Brief)) %>% 
  mutate(Brief = gsub(" serogroup.*" , "",Brief)) %>%  
  mutate(Brief = gsub("*endosymbiont of " , "",Brief)) %>% 
  mutate(Brief = gsub("^ " , "",Brief)) %>% 
  mutate(Brief = gsub("sp.*", "sp.", Brief)) %>% 
  mutate(Brief_cap = toupper(Brief))  %>% 
  select(1,2, 4,Brief_cap)

# Fuzzy match taxonomy file with KEGG
taxo1 <- taxo %>% left_join(., kegglist, by = "Species_cap") %>% 
  group_by(across(c(-ID,-code))) %>% slice_head(n = 1) %>% ungroup() %>% unique()

# Fuzzy match taxonomy file with Uniprot (taxonomy sp. is included in Uniprot sp.)
taxo2 <- taxo1 %>% left_join(., proteome, by = c("Species_cap" = "Brief_cap")) %>% 
  group_by(across(c(-Proteome.Id, -Organism.Id, -Type))) %>% slice_head(n = 1) %>% ungroup() %>% unique()

# Rearrange the taxonomy dataframe 
taxo3 <- taxo2 %>% select(user_genome, Proteome.Id, Organism.Id, code, ID, Species, Domain, Phylum, Class, Order, Family, classification, Type) %>% 
  rename(KEGG.entry = ID, KEGG.org.code = code) %>% 
  mutate(Proteome.Id = if_else(is.na(Proteome.Id), "", Proteome.Id),
         Organism.Id = if_else(is.na(Organism.Id), "", as.character(Organism.Id)),
         KEGG.org.code = if_else(is.na(KEGG.org.code), "", KEGG.org.code),
         KEGG.entry = if_else(is.na(KEGG.entry), "", KEGG.entry)) %>% unique()
write.table(taxo3, "taxonomy.txt",quote = F,sep = "\t",row.names = F,col.names = T)
