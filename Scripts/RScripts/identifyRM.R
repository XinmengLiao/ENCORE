library(dplyr)

sys.args <- commandArgs(trailingOnly = T)

network_file <- sys.args[1]
daa_file <- sys.args[2]
fva_file <- sys.args[3]
output_dir <- sys.args[4]

network <- read.csv(network_file, stringsAsFactors = F,header = T,sep = "\t")
daa <- read.csv(daa_file,header = T, sep = "\t")
fva <- read.csv(fva_file, header = T,sep = "\t")
names(daa)[1:4] <- c("classification","logFC","PValue","FDR")


# merge the daa results with the network 
network1 <- network %>% 
  left_join(., daa %>% select(Vertex = classification, PValue, logFC), by = c("Vertex.A" = "Vertex")) %>% 
  rename(pvalue = PValue, logFC_A = logFC) %>% 
  left_join(., daa %>% select(Vertex = classification, PValue, logFC), by = c("Vertex.B" = "Vertex")) %>% 
  rename(pvalue_B = PValue, logFC_B = logFC) %>% 
  mutate(pvalue = ifelse(is.na(pvalue), pvalue_B, pvalue),
         logFC_A = ifelse(is.na(logFC_A), logFC_B, logFC_A)) %>% 
  select(-logFC_B, -pvalue_B) %>% rename(logFC = logFC_A, PValue= pvalue)


rma <- function(network,pvalue,logFC, reaction){
  
    network <- network %>% 
      filter(!is.na(.data[[pvalue]])) %>% 
      mutate(zscore = qnorm(1 - .data[[pvalue]])) %>%  
      mutate(zscore = ifelse(is.na(zscore), 0, zscore)) %>% 
      mutate(zscore = case_when(zscore == -Inf ~ -15,
                                zscore == Inf ~ 15, TRUE ~ zscore)) 
    
    # For each metabolite, find out the associated microbes
    mets <- unique(c(network %>% filter(Type.A == "Metabolite") %>% pull(Vertex.A),
                     network %>% filter(Type.B == "Metabolite") %>% pull(Vertex.B)))
    
    # calculate weighted zscores
    n_mets <- length(unique(mets))
    n_microbe <- rep(NA, n_mets)
    metZScores <- rep(NA, n_mets)
    metNGenes <- rep(NA, n_mets)
    meanZ <- rep(NA, n_mets)
    stdZ <- rep(NA, n_mets)
    microbe.zscore.list <- list()
    
    for (i in 1:n_mets){
      met.network <- network %>% 
        filter(Vertex.A == mets[i] | Vertex.B == mets[i]) %>% unique()
      metZScores[i] <- sum(met.network$zscore) / sqrt(nrow(met.network))
      meanZ[i] <- mean(met.network$zscore)
      stdZ[i] <- sd(met.network$zscore)
      metNGenes[i] <- mets[i]
      n_microbe[i] <- nrow(met.network)
      microbe.zscore.list[[i]] <- met.network$zscore %>% unlist()
    }
    
    # Filter out NA scores
    valid_mets <- which(!is.na(metZScores))
    mets <- mets[valid_mets]
    n_microbe <- n_microbe[valid_mets]
    metZScores <- metZScores[valid_mets]
    metNGenes <- metNGenes[valid_mets]
    meanZ <- meanZ[valid_mets]
    stdZ <- stdZ[valid_mets]
    
    # Background correction
    set.seed(42)
    n_sim <- 1000
    all_zscores <- unlist(microbe.zscore.list)
    
    # Precompute background stats for each k
    unique_k <- unique(n_microbe)
    background_means <- numeric(length(unique_k))
    background_sds <- numeric(length(unique_k))
    names(background_means) <- names(background_sds) <- unique_k
    
    background_scores_list <- list()
    
    for (j in seq_along(unique_k)) {
      k <- unique_k[j]
      sampled_scores <- replicate(n_sim, {
        sampled_zs <- sample(all_zscores, k, replace = TRUE)
        random_dirs <- ifelse(runif(k) < 0.5, 1, -1)
        directional_zs <- sampled_zs * random_dirs
        sum(directional_zs) / sqrt(k)
      })
      background_means[j] <- mean(sampled_scores)
      background_sds[j] <- sd(sampled_scores)
      background_scores_list[[as.character(k)]] <- sampled_scores
    }
    
    # Apply background correction using precomputed means and sds
    Z_corrected <- numeric(length(metZScores))
    p_corrected <- numeric(length(metZScores))
    
    for (m in seq_along(metZScores)) {
      k <- as.character(n_microbe[m])
      mu_k <- background_means[k]
      sigma_k <- background_sds[k]
      Z_corrected[m] <- (metZScores[m] - mu_k) / sigma_k
      #p_corrected[m] <- 1 - pnorm(abs(Z_corrected[m]))
      p_corrected[m] <- 2 * (1 - pnorm(abs(Z_corrected[m])))
    }
    
    metabolite_scores <- tibble(metabolite = metNGenes,
                                k = n_microbe,
                                weighted_Zdirection = metZScores,
                                weighted_Zdirection_corrected = Z_corrected,
                                p_corrected = p_corrected) %>% 
      as.data.frame() %>% 
      arrange(p_corrected) %>%
      mutate(p_corrected = as.numeric(p_corrected)) %>% 
      left_join(., reaction %>% 
                  select(exchangeRxns, equations) %>% unique(), by = c("metabolite" = "exchangeRxns")) %>% 
      unique() %>%
      mutate(equations = sapply(strsplit(split = " <=> ",equations),`[`,1)) %>% 
      mutate(equations = sapply(strsplit(split = " => ",equations),`[`,1)) %>% 
      mutate(equations = gsub("\\[e\\]$","",equations))
    return(metabolite_scores)
}
  

rm.all <- rma(network1,pvalue = "PValue",logFC = "logFC", reaction = fva)

network.increase <- network1 %>% filter(logFC > 0)
rm.up <- rma(network.increase, pvalue = "PValue", logFC = "logFC", reaction = fva)

network.decrease <- network1 %>% filter(logFC < 0) 
rm.down <- rma(network.decrease, pvalue = "PValue", logFC = "logFC", reaction = fva)

write.table(rm.all, file = file.path(output_dir,"reporter_metabolites.txt"), sep = "\t", row.names = F, quote = F)
write.table(rm.up, file = file.path(output_dir,"reporter_metabolites_increased.txt"), sep = "\t", row.names = F, quote = F)
write.table(rm.down, file = file.path(output_dir,"reporter_metabolites_decreased.txt"), sep = "\t", row.names = F, quote = F)