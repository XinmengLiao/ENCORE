library(dplyr)
library(stringr)
library(webchem)

#metname <- read.csv("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/Feces_10028_Visit1_S98/bin.1.s/data/bin1s.unique_metabolites.txt",header = F,sep = "\t")
samples <- list.files("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/")
#done <- c("Feces_10010_Visit1_S106","Feces_20004_Visit1_S5","Feces_20004_Visit3_S9","Feces_20030_Visit1_S8","Feces_20030_Visit3_S85","Feces_20054_Visit3_S14","Feces_20073_Visit3_S2")
#samples <- setdiff(samples,done)
sample <- unlist(df)

name.edit <- function(data){
  metname = data
  chemical_formula_pattern <- "(\\sC[0-9]+H[0-9]+[A-Za-z0-9]*)$"
  #chemical_formula_pattern <- "(\s[A-Z][a-z]?[0-9]*([A-Z][a-z]?[0-9]*)*)$"
  metname <- metname %>%
    mutate(
      has_chemical_formula = str_detect(col, chemical_formula_pattern)) %>% 
    mutate(col = ifelse(
      has_chemical_formula, 
      str_remove(col, chemical_formula_pattern),
      col)) %>% 
    mutate(col = if_else(has_chemical_formula, gsub("^ R ((.*?)[A-Za-z0-9]+)$", "(R)\\1", col), col),
           col = if_else(has_chemical_formula, gsub("^ S ((.*?)[A-Za-z0-9]+)$", "(S)\\1", col), col),
           col = if_else(has_chemical_formula, gsub("(\\d) (\\d)", "\\1,\\2", col), col)) %>% 
    mutate(col = if_else(has_chemical_formula, gsub(" ", "-", col), col)) %>% 
    mutate(col = if_else(!has_chemical_formula, gsub("^D ", "D-", col),col),
           col = if_else(!has_chemical_formula, gsub("^L ", "L-", col),col),
           col = if_else(!has_chemical_formula, gsub("^N ", "N-", col),col),
           col = gsub("^N N", "N,N",col),
           col = gsub("^ (\\d) ","\\1-",col),
           col = gsub("^ ([RS0-9]+) ([RS0-9]+) ", "(\\1,\\2)", col),
           col = gsub("\\b([A-Za-z0-9]+) \\1\\b", "\\1", col),
           col = if_else(has_chemical_formula, gsub("  ([a-zA-Z0-9]+) (.*)([CHONPS]|Fe|[0-9])+$", " (\\1)\\2\\3", col),col)) %>% 
    select(-has_chemical_formula) 
  return(metname)
}
pubchem.api <- function(metname.rest){
  # Retrieve CID from PubChem via API
  all_cids <- data.frame()  
  n <- length(metname.rest$col)
  chunk_size <- 50
  for (cidnum in seq(1, n, by = chunk_size)) {  
    print(cidnum)
    print("Retrieving CID")
    end <- min(cidnum + chunk_size - 1, n)
    chunk <- metname.rest$col[cidnum:end]
    cids <- get_cid(chunk)
    all_cids <- rbind(all_cids, cids)
  }
  all_cids <- unique(all_cids)
  
  
  index <- unique(all_cids$cid)
  if(length(index) > 1 & any(!is.na(index))){
      # Retrieve SMILEs from PubChem via API
      all_smiles <- data.frame() 
      n <- nrow(all_cids)
      chunk_size <- 100
      for (smilenum in seq(1, n, by = chunk_size)) {
        print("Retrieving SMILEs")
        print(smilenum)
        end <- min(smilenum + chunk_size - 1, n)
        chunk <- all_cids$cid[smilenum:end]
        smiles <- pc_prop(chunk, properties = "CanonicalSMILES") 
        all_smiles <- rbind(all_smiles, smiles)
      }
      all_smiles <- unique(all_smiles)
      
      # Merged
      all_smiles$CID <- as.character(all_smiles$CID)
      cid.smiles <- metname.rest %>% 
        select(1:2) %>% 
        left_join(., all_cids, by = c("col" = "query")) %>% 
        left_join(., all_smiles, by = c("cid" = "CID")) %>% unique() %>% 
        mutate(source = "metName.DB")
      cid.smiles <- cid.smiles[!duplicated(cid.smiles$Org.name),]
      return(cid.smiles)
  }else{
    cid.smiles <- metname.rest %>% 
      mutate(source = "metName.DB")
    cid.smiles <- cid.smiles[!duplicated(cid.smiles$Org.name),]
    return(cid.smiles)
    }
}


for (j in samples){
  print(j)
  bins <- list.files(paste0("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/",j)) 
  bins <- bins[grep("bin",bins)]
  
  for (i in bins){
    print(c(j,i))
    # Load metName database
    met.inhouse.db.path <- "/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/metName.inhouseDB.txt"
    metname.inhouse.database <- read.csv(met.inhouse.db.path,header = T,sep = "\t")
    
    binID = gsub("\\.","",i)
    metname <- read.csv(paste0("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/",j,"/",
                               i,"/data/",binID,".unique_metabolites.txt"),header = F,sep = "\t")
    colnames(metname) = "Org.name"
    metname$col = metname$Org.name
    metname <- name.edit(data = metname)
    
    # Compare with the MetSmiles DB
    metname <- metname %>% 
      left_join(., metname.inhouse.database, by = c("Org.name","col"))
    
    metname.exist <- metname %>% filter(!is.na(source))
    metname.rest <- metname %>% filter(is.na(source))
    
    if(nrow(metname.rest) > 0){
      metname.rest <- pubchem.api(metname.rest)
      metName.all <- rbind(metname.exist,metname.rest) %>% unique() 
      write.table(metName.all,paste0("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/",j,"/",
                                     i,"/data/",binID,".unique_metabolites.smiles.txt"),quote = F,sep = "\t",row.names = F)
      metName.all.final <- metName.all %>% 
        select(Org.name,CanonicalSMILES) %>% 
        mutate(CanonicalSMILES = if_else(is.na(CanonicalSMILES),"",CanonicalSMILES))
      write.table(metName.all.final,
                  paste0("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/",j,"/",
                                                                          i,"/data/smilesDB.tsv"),quote = F,sep = "\t",
                  row.names = F,col.names = F)
      
      # Accmulate to database
      metname.inhouse.database <- rbind(metname.inhouse.database, metname.rest) %>% unique()
      write.table(metname.inhouse.database, met.inhouse.db.path,quote = F,sep = "\t",row.names = F)
    }else{
      write.table(metname.exist,paste0("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/",j,"/",
                                       i,"/data/",binID,".unique_metabolites.smiles.txt"),quote = F,sep = "\t",row.names = F)
      metName.all.final <- metname.exist %>% 
        select(Org.name,CanonicalSMILES) %>% 
        mutate(CanonicalSMILES = if_else(is.na(CanonicalSMILES),"",CanonicalSMILES))
      write.table(metName.all.final %>% select(Org.name,CanonicalSMILES),
                  paste0("/Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/",j,"/",
                         i,"/data/smilesDB.tsv"),quote = F,sep = "\t",
                  row.names = F,col.names = F)
    }
  }
}





