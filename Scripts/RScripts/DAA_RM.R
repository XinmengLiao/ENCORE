library(dplyr)
library(tidyr)
library(Maaslin2)

# Differential abundance analysis (DAA) of species with MaAsLin2, microbe-metabolite network
# from the mandatory exchange reactions, and directional reporter metabolites (PredCMB-style)
# Usage:
#   Rscript DAA_RM.R <abundance_file> <taxo_dir> <metadata_file> <fva_file> <output_dir> \
#                    [fixed_effect] [random_effect]
#
# abundance_file : tab-separated, columns Sample, Bin, Mapped_reads, Bin_length
# taxo_dir       : folder with <sample>_gtdbtk_summary.tsv files from GTDB-Tk
# metadata_file  : tab-separated, first column is the sample ID (same as Sample in abundance_file)
# fva_file       : mandatory exchange reactions from ENCORE.sh --fva (fva/fva_mandatory_exchange.txt),
#                  columns exchangeRxns, equations, LB_exc, UB_exc, Direction, Sample, Bin
# output_dir     : output folder
# fixed_effect   : metadata column compared between two groups (default: Group); more columns can be
#                  added as covariates, comma separated, the first one is used for the reporter metabolites
# random_effect  : metadata column(s) for repeated measures, e.g. subject for paired samples (default: NA)
#
# Output in output_dir:
#   maaslin2/                         MaAsLin2 output
#   taxoID.txt                        species ID used in MaAsLin2 and its GTDB-Tk classification
#   daa.txt                           DAA results (classification, logFC, PValue, FDR)
#   MMNetwork.txt                     microbe-metabolite network
#   reporter_metabolites.txt          reporter metabolites in the entire network
#   reporter_metabolites_uptake.txt   reporter metabolites in the uptake subnetwork
#   reporter_edge_scores.txt          signed Z-score of each microbe-metabolite edge (entire network)

sys.args <- commandArgs(trailingOnly = TRUE)
abundance_file <- sys.args[1]
taxo_dir <- sys.args[2]
metadata_file <- sys.args[3]
fva_file <- sys.args[4]
output_dir <- sys.args[5]
fixed_effect <- ifelse(is.na(sys.args[6]), "Group", sys.args[6])
random_effect <- ifelse(is.na(sys.args[7]), "NA", sys.args[7])

fixed_effect <- unlist(strsplit(fixed_effect, ","))
random_effect <- if (random_effect == "NA") NULL else unlist(strsplit(random_effect, ","))
dir.create(output_dir, showWarnings = F, recursive = T)

####-------------------------------------------####
#### 1. DAA with MaAsLin2
####-------------------------------------------####

# read taxonomy files from all samples
taxo.all <- data.frame()
taxo_files <- list.files(taxo_dir, pattern = "_gtdbtk_summary")
print("Now reading taxonomy files from all samples. ")
for (i in taxo_files){
  sampleid = unlist(strsplit(i,"_gtdbtk_summary"))[1]
  taxo.tmp <- read.csv(file.path(taxo_dir,i),header = T,sep = "\t") %>%
    mutate(Sample = sampleid) %>% select(Sample, user_genome, classification) %>% unique()
  taxo.all <- rbind(taxo.all, taxo.tmp)
}
taxo.all <- taxo.all %>%
  mutate(user_genome = gsub("\\.","",user_genome)) %>%
  filter(grepl(";s__.+", classification))

# read abundance of all reassembled bins
map.summary <- read.csv(abundance_file, header = T, sep = "\t") %>%
  mutate(user_genome = gsub("\\.","",Bin))

metadata <- read.csv(metadata_file, header = T, sep = "\t", check.names = F)
colnames(metadata)[1] <- "Sample"
map.summary <- map.summary %>% filter(Sample %in% metadata$Sample)
print(paste0("Number of samples used for DAA: ", length(unique(map.summary$Sample))))

# calculate the relative abundance, and sum bins of the same species
abund <- map.summary %>%
  group_by(Sample) %>%
  mutate(Relative_abundance = (Mapped_reads / Bin_length) / sum(Mapped_reads / Bin_length, na.rm = TRUE)) %>%
  ungroup() %>%
  inner_join(taxo.all, by = c("Sample", "user_genome")) %>%
  group_by(Sample, classification) %>%
  summarise(Relative_abundance = sum(Relative_abundance), .groups = "drop")

# use short species IDs as feature names, since MaAsLin2 changes special characters in the names
taxo.id <- data.frame(classification = sort(unique(abund$classification))) %>%
  mutate(SpeciesID = paste0("Species", row_number()))
write.table(taxo.id, file.path(output_dir, "taxoID.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

feature_table <- abund %>%
  left_join(taxo.id, by = "classification") %>%
  select(Sample, SpeciesID, Relative_abundance) %>%
  pivot_wider(names_from = SpeciesID, values_from = Relative_abundance, values_fill = 0) %>%
  as.data.frame()
rownames(feature_table) <- feature_table$Sample
feature_table$Sample <- NULL

metadata <- metadata %>% filter(Sample %in% rownames(feature_table)) %>% as.data.frame()
rownames(metadata) <- metadata$Sample
metadata$Sample <- NULL
for (i in fixed_effect) metadata[[i]] <- as.factor(metadata[[i]])
for (i in random_effect) metadata[[i]] <- as.character(metadata[[i]])
if (nlevels(metadata[[fixed_effect[1]]]) != 2){
  stop(fixed_effect[1], " must have two groups for the reporter metabolites, found: ",
       paste(levels(metadata[[fixed_effect[1]]]), collapse = ", "))
}

fit_data <- Maaslin2(
  input_data = feature_table,
  input_metadata = metadata,
  output = file.path(output_dir, "maaslin2"),
  fixed_effects = fixed_effect,
  random_effects = random_effect, # paired
  min_prevalence = 0.05,
  min_abundance = 0
)

# columns classification, logFC, PValue, FDR; only the comparison of the first fixed effect
masslin.res <- read.csv(file.path(output_dir, "maaslin2", "all_results.tsv"), sep = "\t", header = T) %>%
  filter(metadata == fixed_effect[1]) %>%
  rename(SpeciesID = feature, logFC = coef, PValue = pval, FDR = qval) %>%
  left_join(taxo.id, by = "SpeciesID") %>%
  select(classification, logFC, PValue, FDR, SpeciesID, metadata, value)
write.table(masslin.res, file.path(output_dir, "daa.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
print(paste0("Comparison: ", fixed_effect[1], " ", unique(masslin.res$value), " vs ", levels(metadata[[fixed_effect[1]]])[1]))

####-------------------------------------------####
#### 2. Microbe-metabolite network
####-------------------------------------------####

# mandatory exchange reactions (fva_data_manage in ENCORE.sh): UB_exc < 0 Uptake, LB_exc > 0 Secretion
fva.all <- read.csv(fva_file, header = T, sep = "\t") %>%
  mutate(user_genome = gsub("\\.","",Bin),
         direction = if_else(Direction == "Uptake", "Uptake", "Secrete")) %>%
  unique()
print(paste0("Number of samples with FVA results: ", length(unique(fva.all$Sample))))

build.fva.network <- function(daa_result, fva.org, pvalue, padj, log2fc) {

  fva.annotated <- fva.org %>%
    left_join(taxo.all, by = c("Sample", "user_genome")) %>%
    filter(!is.na(classification)) %>%
    left_join(daa_result, by = "classification") %>%
    filter(!is.na(.data[[log2fc]])) %>%
    distinct()

  secrete <- fva.annotated %>% filter(direction == "Secrete") %>%
    select(classification, exchangeRxns, direction, all_of(c(pvalue, padj, log2fc))) %>%
    unique() %>%
    rename(Vertex.A = classification, Vertex.B = exchangeRxns) %>%
    mutate(Type.A = "Microbe", Type.B = "Metabolite")
  network <- fva.annotated %>% filter(direction != "Secrete") %>%
    select(exchangeRxns, classification, direction, all_of(c(pvalue, padj, log2fc))) %>%
    unique() %>%
    rename(Vertex.A = exchangeRxns, Vertex.B = classification) %>%
    mutate(Type.A = "Metabolite", Type.B = "Microbe") %>%
    rbind(., secrete) %>%
    unique() %>%
    filter(!is.na(.data[[padj]])) %>%
    filter(Vertex.A != Vertex.B) %>%
    # a species that both takes up and secretes a metabolite (in different bins) is removed for it
    group_by(edge_group = paste(pmin(Vertex.A, Vertex.B), pmax(Vertex.A, Vertex.B))) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    select(-edge_group) %>%
    select(Vertex.A, Type.A, Vertex.B, Type.B, everything())

  microbe.num <- network %>%
    mutate(microbe = if_else(Type.B == "Microbe", Vertex.B, Vertex.A)) %>%
    pull(microbe) %>% unique()
  print(paste0("Number of microbes in the network: ", length(microbe.num)))

  return(network)
}

network <- build.fva.network(daa_result = masslin.res %>% select(classification, logFC, PValue, FDR),
                             fva.org = fva.all, pvalue = "PValue", padj = "FDR", log2fc = "logFC")
network.uptake <- network %>% filter(direction == "Uptake")
write.table(network, file.path(output_dir, "MMNetwork.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

####-------------------------------------------####
#### 3. Directional reporter metabolites 
####-------------------------------------------####
# PredCMB algorithm with microbial species in place of enzymatic gene families. No P-value or logFC
# threshold, no abundance weighting. Steps:
#   1. Two-sided species P -> abundance-directional P: increased = P/2, decreased = 1 - P/2 (logFC = 0: Z = 0)
#   2. Directional P -> Z = qnorm(1 - P_directional); only infinite Z are replaced by +/-max_abs_z (7)
#   3. Exchange direction: Secrete = +1, Uptake = -1; edge_Z = z_abundance * reaction_sign
#   4. Metabolite score = sum(edge_Z) / k (mean), k = number of unique species linked to the metabolite
#   5. For each metabolite, sample k unique background species n_sim times (no replacement),
#      keeping its numbers of secretion and uptake edges
#   6. Standardize the score with the null mean/sd, P from the sign-specific normal tail (p_corrected)
#   7. BH-FDR across metabolites; two-sided empirical permutation P/FDR reported as sensitivity result
rma.predcmb <- function(network, pvalue, logFC, reaction,
                        n_sim = 1000, seed = 42, max_abs_z = 7) {

  # microbe and metabolite of each edge
  target_edges <- network %>%
    mutate(
      microbe = case_when(Type.A == "Microbe" ~ Vertex.A, Type.B == "Microbe" ~ Vertex.B, TRUE ~ NA_character_),
      metabolite = case_when(Type.A == "Metabolite" ~ Vertex.A, Type.B == "Metabolite" ~ Vertex.B, TRUE ~ NA_character_)
    ) %>%
    filter(!is.na(microbe), !is.na(metabolite), direction %in% c("Secrete", "Uptake"))

  if (anyNA(target_edges[[pvalue]]) || anyNA(target_edges[[logFC]])) {
    stop("Missing values found in ", pvalue, " or ", logFC, ".")
  }
  if (any(!is.finite(target_edges[[pvalue]])) ||
      any(target_edges[[pvalue]] < 0 | target_edges[[pvalue]] > 1)) {
    stop(pvalue, " must contain finite values between 0 and 1.")
  }
  if (any(!is.finite(target_edges[[logFC]]))) {
    stop(logFC, " must contain finite values.")
  }

  # Step 1-3
  target_edges <- target_edges %>%
    mutate(
      p_input = .data[[pvalue]],
      abundance_change = case_when(.data[[logFC]] > 0 ~ "Increased", .data[[logFC]] < 0 ~ "Decreased", TRUE ~ "Unchanged"),
      p_directional_abundance = case_when(.data[[logFC]] > 0 ~ p_input / 2, .data[[logFC]] < 0 ~ 1 - p_input / 2, TRUE ~ 0.5),
      z_abundance_uncapped = qnorm(1 - p_directional_abundance),
      z_abundance = case_when(z_abundance_uncapped == Inf ~ max_abs_z,
                              z_abundance_uncapped == -Inf ~ -max_abs_z,
                              TRUE ~ z_abundance_uncapped),
      reaction_sign = if_else(direction == "Secrete", 1, -1),
      edge_Z = z_abundance * reaction_sign
    ) %>%
    distinct(metabolite, microbe, direction, .keep_all = TRUE)

  # the same microbe-metabolite pair cannot be both uptake and secretion
  direction_conflicts <- target_edges %>%
    count(metabolite, microbe, name = "n_directions") %>%
    filter(n_directions > 1)
  if (nrow(direction_conflicts) > 0) {
    stop("Some microbe-metabolite pairs have both Secrete and Uptake directions.")
  }
  target_edges <- target_edges %>% distinct(metabolite, microbe, .keep_all = TRUE)

  # background of unique species, no weights
  microbial_background <- target_edges %>%
    select(microbe, p_input, p_directional_abundance, z_abundance, all_of(logFC), all_of(pvalue)) %>%
    distinct(microbe, .keep_all = TRUE)

  metabolite_edges <- target_edges %>%
    group_by(metabolite) %>%
    group_split(.keep = TRUE)

  # Step 4-6
  set.seed(seed)
  score_one <- function(edges) {
    k <- nrow(edges)
    n_secrete <- sum(edges$reaction_sign == 1)
    n_uptake <- sum(edges$reaction_sign == -1)
    direction_template <- c(rep(1, n_secrete), rep(-1, n_uptake))
    observed_z <- mean(edges$edge_Z)

    if (k > nrow(microbial_background)) {
      return(tibble(
        metabolite = edges$metabolite[1], k = k, n_secrete = n_secrete,
        n_uptake = n_uptake, directional_Z_raw = observed_z,
        null_mean = NA_real_, null_sd = NA_real_, directional_Z_null_adjusted = NA_real_,
        p_corrected = NA_real_, p_empirical = NA_real_
      ))
    }

    null_scores <- replicate(n_sim, {
      sampled <- microbial_background %>% slice_sample(n = k, replace = FALSE)
      mean(sampled$z_abundance * direction_template)
    })

    null_mean <- mean(null_scores)
    null_sd <- sd(null_scores)
    z_adjusted <- if (is.finite(null_sd) && null_sd > 0) (observed_z - null_mean) / null_sd else NA_real_
    p_corrected <- if (is.na(z_adjusted)) {
      NA_real_
    } else if (z_adjusted > 0) {
      1 - pnorm(z_adjusted)
    } else {
      pnorm(z_adjusted)
    }
    p_lower <- (1 + sum(null_scores <= observed_z)) / (n_sim + 1)
    p_upper <- (1 + sum(null_scores >= observed_z)) / (n_sim + 1)
    p_empirical <- min(1, 2 * min(p_lower, p_upper))
    tibble(
      metabolite = edges$metabolite[1], k = k,
      n_secrete = n_secrete, n_uptake = n_uptake,
      directional_Z_raw = observed_z,
      null_mean = null_mean, null_sd = null_sd,
      directional_Z_null_adjusted = z_adjusted,
      p_corrected = p_corrected,
      p_empirical = p_empirical
    )
  }

  # Step 7
  metabolite_scores <- bind_rows(lapply(metabolite_edges, score_one)) %>%
    mutate(
      predicted_direction = case_when(directional_Z_raw > 0 ~ "Increased",
                                      directional_Z_raw < 0 ~ "Decreased",
                                      TRUE ~ "Unchanged"),
      FDR_corrected = p.adjust(p_corrected, method = "BH"),
      FDR_empirical = p.adjust(p_empirical, method = "BH")
    ) %>%
    left_join(
      reaction %>% select(exchangeRxns, equations) %>%
        filter(!is.na(equations)) %>% distinct(exchangeRxns, .keep_all = TRUE),
      by = c("metabolite" = "exchangeRxns")
    ) %>%
    mutate(
      equations = sub(" <=> .*", "", equations),
      equations = sub(" => .*", "", equations),
      equations = sub("\\[e\\]$", "", equations)
    ) %>%
    arrange(FDR_corrected, p_corrected)

  edge_scores <- target_edges %>%
    mutate(predicted_edge_direction = case_when(edge_Z > 0 ~ "Increased",
                                                edge_Z < 0 ~ "Decreased",
                                                TRUE ~ "Unchanged"))

  list(all.met = metabolite_scores, edge.scores = edge_scores)
}

print("Now identifying reporter metabolites in the entire network. ")
rm.all <- rma.predcmb(network, pvalue = "PValue", logFC = "logFC", reaction = fva.all)
print("Now identifying reporter metabolites in the uptake subnetwork. ")
rm.uptake <- rma.predcmb(network.uptake, pvalue = "PValue", logFC = "logFC", reaction = fva.all)

write.table(rm.all$all.met, file.path(output_dir, "reporter_metabolites.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(rm.uptake$all.met, file.path(output_dir, "reporter_metabolites_uptake.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(rm.all$edge.scores, file.path(output_dir, "reporter_edge_scores.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
print(paste0("Reporter metabolites (p_corrected < 0.05): ", sum(rm.all$all.met$p_corrected < 0.05, na.rm = T),
             " in the entire network, ", sum(rm.uptake$all.met$p_corrected < 0.05, na.rm = T), " in the uptake subnetwork"))
print(paste0("Results are stored in ", output_dir))
