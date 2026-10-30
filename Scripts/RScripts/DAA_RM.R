library(dplyr)
library(tidyr)
library(Maaslin2)

# Differential abundance analysis (DAA) of species with MaAsLin2, microbe-metabolite network
# from the mandatory exchange reactions, and directional reporter metabolites (PredCMB-style,
# weighted by the mean relative abundance of each species)
# Usage:
#   Rscript DAA_RM.R <abundance_file> <taxo_dir> <metadata_file> <fva_file> <output_dir> \
#                    [fixed_effect] [checkm_file] [completeness] [contamination]
#
# abundance_file : tab-separated, columns Sample, Bin, Mapped_reads, Bin_length; one row per bin of
#                  all samples (ENCORE.sh --reporter collects it from the .map files of --abundance)
# taxo_dir       : folder with <sample>_gtdbtk_summary.tsv files from GTDB-Tk
# metadata_file  : tab-separated, first column is the sample ID (same as Sample in abundance_file)
# fva_file       : mandatory exchange reactions from ENCORE.sh --fva (fva/fva_mandatory_exchange.txt),
#                  columns exchangeRxns, equations, LB_exc, UB_exc, Direction, Sample, Bin
# output_dir     : output folder
# fixed_effect   : metadata column with the two groups compared (default: Group); MaAsLin2 is run with
#                  this column only, without covariates or random effects (no paired design)
# checkm_file    : tab-separated, columns Sample, Bin, completeness, contamination (CheckM of the reassembled
#                  bins, collected by ENCORE.sh --reporter); default: NA, all bins are used
# completeness   : only bins with completeness > this value (default: 70)
# contamination  : and contamination < this value (default: 10)
#
# Output in output_dir:
#   maaslin2/                         MaAsLin2 output
#   taxoID.txt                        species ID used in MaAsLin2 and its GTDB-Tk classification
#   daa.txt                           DAA results (classification, logFC, PValue, FDR)
#   MMNetwork.txt                     microbe-metabolite network
#   species_abundance.txt             relative abundance of each species in each sample (bins summed)
#   reporter_metabolites.txt          reporter metabolites in the entire network
#   reporter_metabolites_uptake.txt   reporter metabolites in the uptake subnetwork
#   reporter_edge_scores.txt          signed Z-score and abundance weight of each microbe-metabolite edge
#                                     (entire network)
#   reporter_microbial_background.txt background species and their abundance weights (both networks)

sys.args <- commandArgs(trailingOnly = TRUE)
abundance_file <- sys.args[1]
taxo_dir <- sys.args[2]
metadata_file <- sys.args[3]
fva_file <- sys.args[4]
output_dir <- sys.args[5]
fixed_effect <- ifelse(is.na(sys.args[6]), "Group", sys.args[6])
checkm_file <- ifelse(is.na(sys.args[7]), "NA", sys.args[7])
completeness_cutoff <- ifelse(is.na(sys.args[8]), 70, as.numeric(sys.args[8]))
contamination_cutoff <- ifelse(is.na(sys.args[9]), 10, as.numeric(sys.args[9]))

if (grepl(",", fixed_effect)) stop("fixed_effect must be a single metadata column, found: ", fixed_effect)
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

# keep the bins that pass the CheckM cutoffs; the relative abundance is calculated among these bins
if (checkm_file != "NA") {
  checkm <- read.csv(checkm_file, header = T, sep = "\t") %>%
    filter(completeness > completeness_cutoff, contamination < contamination_cutoff)
  map.summary <- map.summary %>% semi_join(checkm, by = c("Sample", "Bin"))
  print(paste0("Bins with completeness > ", completeness_cutoff, " and contamination < ", contamination_cutoff,
               ": ", nrow(map.summary)))
}
print(paste0("Number of samples used for DAA: ", length(unique(map.summary$Sample))))

# calculate the relative abundance, and sum bins of the same species
abund <- map.summary %>%
  group_by(Sample) %>%
  mutate(Relative_abundance = (Mapped_reads / Bin_length) / sum(Mapped_reads / Bin_length, na.rm = TRUE)) %>%
  ungroup() %>%
  inner_join(taxo.all, by = c("Sample", "user_genome")) %>%
  group_by(Sample, classification) %>%
  summarise(Relative_abundance = sum(Relative_abundance), .groups = "drop")
write.table(abund, file.path(output_dir, "species_abundance.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

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
if (!fixed_effect %in% colnames(metadata)) stop(fixed_effect, " is not a column of ", metadata_file)
metadata <- metadata[, fixed_effect, drop = FALSE]
metadata[[fixed_effect]] <- as.factor(metadata[[fixed_effect]])
if (nlevels(metadata[[fixed_effect]]) != 2){
  stop(fixed_effect, " must have two groups for the reporter metabolites, found: ",
       paste(levels(metadata[[fixed_effect]]), collapse = ", "))
}

fit_data <- Maaslin2(
  input_data = feature_table,
  input_metadata = metadata,
  output = file.path(output_dir, "maaslin2"),
  fixed_effects = fixed_effect, # group only: no covariates, no random effects
  min_prevalence = 0.05,
  min_abundance = 0
)

# columns classification, logFC, PValue, FDR
masslin.res <- read.csv(file.path(output_dir, "maaslin2", "all_results.tsv"), sep = "\t", header = T) %>%
  filter(metadata == fixed_effect) %>%
  rename(SpeciesID = feature, logFC = coef, PValue = pval, FDR = qval) %>%
  left_join(taxo.id, by = "SpeciesID") %>%
  select(classification, logFC, PValue, FDR, SpeciesID, metadata, value)
write.table(masslin.res, file.path(output_dir, "daa.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
print(paste0("Comparison: ", fixed_effect, " ", unique(masslin.res$value), " vs ", levels(metadata[[fixed_effect]])[1]))

####-------------------------------------------####
#### 2. Microbe-metabolite network
####-------------------------------------------####

# mandatory exchange reactions (fva_data_manage in ENCORE.sh): UB_exc < 0 Uptake, LB_exc > 0 Secretion
fva.all <- read.csv(fva_file, header = T, sep = "\t") %>%
  mutate(user_genome = gsub("\\.","",Bin),
         direction = if_else(Direction == "Uptake", "Uptake", "Secrete")) %>%
  unique()
# use the same CheckM-passing bins as the abundance, so every edge comes from a bin counted in the DAA
if (checkm_file != "NA") {
  fva.all <- fva.all %>% semi_join(checkm, by = c("Sample", "Bin"))
  print(paste0("FVA bins with completeness > ", completeness_cutoff, " and contamination < ", contamination_cutoff,
               ": ", nrow(distinct(fva.all, Sample, Bin))))
}
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
# PredCMB algorithm with microbial species in place of enzymatic gene families, weighted by species
# abundance. No P-value or logFC threshold. Steps:
#   1. Two-sided species P -> abundance-directional P: increased = P/2, decreased = 1 - P/2 (logFC = 0: Z = 0)
#   2. Directional P -> Z = qnorm(1 - P_directional); only infinite Z are replaced by +/-max_abs_z (7)
#   3. Exchange direction: Secrete = +1, Uptake = -1; edge_Z = z_abundance * reaction_sign
#   4. Abundance weight of a species: relative abundance summed over its bins within each sample, then
#      averaged over all samples (nondetection = 0), and normalized by the sum over all background species
#   5. Metabolite score = sum(weight * edge_Z) / k, k = number of unique species linked to the metabolite
#   6. For each metabolite, sample k unique background species (Z and weight paired) n_sim times
#      (no replacement), keeping its numbers of secretion and uptake edges
#   7. Standardize the score with the null mean/sd; p_corrected is the right-tailed normal P for a
#      positive adjusted Z and the left-tailed normal P for a negative one (as in PredCMB)
#   8. BH-FDR across all metabolites
# A metabolite linked to all background species (k = N) with only uptake or only secretion edges has
# a null without variance (every draw is the full background with the same directions): its adjusted
# Z-score is NA and p_corrected = 1, and it is kept in the BH correction.
rma.predcmb <- function(network, pvalue, logFC, reaction, abundance_data,
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

  # Step 4: mean relative abundance of each species over all samples, nondetection counted as 0
  if (any(!is.finite(abundance_data$relative_abundance)) ||
      any(abundance_data$relative_abundance < 0 | abundance_data$relative_abundance > 1 + 1e-8)) {
    stop("Species relative abundance must be finite values between 0 and 1.")
  }
  n_abundance_samples <- n_distinct(abundance_data$sample)
  abundance_weights <- abundance_data %>%
    group_by(sample, microbe) %>%
    summarise(sample_abundance = sum(relative_abundance), .groups = "drop") %>%
    group_by(microbe) %>%
    summarise(abundance_weight = sum(sample_abundance) / n_abundance_samples, .groups = "drop")
  missing_microbes <- setdiff(unique(target_edges$microbe), abundance_weights$microbe)
  if (length(missing_microbes) > 0) {
    stop("No abundance for ", length(missing_microbes), " network species:\n",
         paste(missing_microbes, collapse = "\n"))
  }
  target_edges <- target_edges %>% left_join(abundance_weights, by = "microbe")
  if (any(target_edges$abundance_weight <= 0)) {
    stop("Network species must have a positive mean abundance.")
  }

  # background of unique species, weights normalized by the sum over the background
  microbial_background <- target_edges %>%
    select(microbe, abundance_weight, p_input, p_directional_abundance, z_abundance,
           all_of(logFC), all_of(pvalue)) %>%
    distinct(microbe, .keep_all = TRUE)
  global_weight_total <- sum(microbial_background$abundance_weight)
  target_edges$score_weight <- target_edges$abundance_weight / global_weight_total
  microbial_background$score_weight <- microbial_background$abundance_weight / global_weight_total

  metabolite_edges <- target_edges %>%
    group_by(metabolite) %>%
    group_split(.keep = TRUE)

  # Step 5-7
  background_z <- microbial_background$z_abundance
  background_weights <- microbial_background$score_weight
  n_background <- length(background_z)
  set.seed(seed)
  score_one <- function(edges) {
    k <- nrow(edges)
    n_secrete <- sum(edges$reaction_sign == 1)
    n_uptake <- sum(edges$reaction_sign == -1)
    direction_template <- c(rep(1, n_secrete), rep(-1, n_uptake))
    if (k > n_background) {
      stop("Cannot draw ", k, " unique microbes from a background of ", n_background, ".")
    }
    observed_z <- sum(edges$edge_Z * edges$score_weight) / k

    # k background species per draw, Z and weight kept paired
    sampled_idx <- matrix(replicate(n_sim, sample.int(n_background, k, replace = FALSE)), nrow = k)
    sampled_z <- matrix(background_z[sampled_idx], nrow = k)
    sampled_weights <- matrix(background_weights[sampled_idx], nrow = k)
    null_scores <- colSums(sampled_z * sampled_weights * direction_template) / k

    null_mean <- mean(null_scores)
    null_sd <- sd(null_scores)
    z_adjusted <- if (is.finite(null_sd) && null_sd > 0) (observed_z - null_mean) / null_sd else NA_real_
    # PredCMB: upper normal tail for a positive adjusted Z, lower tail for a negative one
    p_corrected <- if (is.na(z_adjusted)) NA_real_ else if (z_adjusted > 0) 1 - pnorm(z_adjusted) else pnorm(z_adjusted)
    # all background species in one direction only: no null variance, adjusted Z undefined, P = 1
    # (the draws above are still made, so the random numbers of the other metabolites do not change)
    if (k == n_background && (n_secrete == 0 || n_uptake == 0)) {
      null_mean <- observed_z
      null_sd <- 0
      z_adjusted <- NA_real_
      p_corrected <- 1
    }
    tibble(
      metabolite = edges$metabolite[1], k = k, n_background = n_background,
      n_secrete = n_secrete, n_uptake = n_uptake,
      directional_Z_raw = observed_z,
      null_mean = null_mean, null_sd = null_sd,
      directional_Z_null_adjusted = z_adjusted,
      p_corrected = p_corrected
    )
  }

  # Step 8: BH over all metabolites, a P that cannot be calculated counted as 1
  metabolite_scores <- bind_rows(lapply(metabolite_edges, score_one)) %>%
    mutate(
      predicted_direction = case_when(directional_Z_raw > 0 ~ "Increased",
                                      directional_Z_raw < 0 ~ "Decreased",
                                      TRUE ~ "Unchanged"),
      FDR_corrected = p.adjust(replace(p_corrected, is.na(p_corrected), 1), method = "BH")
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

  list(all.met = metabolite_scores, edge.scores = edge_scores,
       microbial.background = microbial_background, n_abundance_samples = n_abundance_samples)
}

# relative abundance of each species in every sample used in the DAA (bins of a species summed)
species_abundance <- abund %>%
  filter(Sample %in% rownames(metadata)) %>%
  transmute(sample = Sample, microbe = classification, relative_abundance = Relative_abundance)

print("Now identifying reporter metabolites in the entire network. ")
rm.all <- rma.predcmb(network, pvalue = "PValue", logFC = "logFC", reaction = fva.all,
                      abundance_data = species_abundance)
print("Now identifying reporter metabolites in the uptake subnetwork. ")
rm.uptake <- rma.predcmb(network.uptake, pvalue = "PValue", logFC = "logFC", reaction = fva.all,
                         abundance_data = species_abundance)
print(paste0("Number of samples used for the abundance weights: ", rm.all$n_abundance_samples))

write.table(rm.all$all.met, file.path(output_dir, "reporter_metabolites.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(rm.uptake$all.met, file.path(output_dir, "reporter_metabolites_uptake.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(rm.all$edge.scores, file.path(output_dir, "reporter_edge_scores.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(bind_rows(rm.all$microbial.background %>% mutate(network = "full"),
                      rm.uptake$microbial.background %>% mutate(network = "uptake")),
            file.path(output_dir, "reporter_microbial_background.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

one_direction_all <- function(res) sum(res$k == res$n_background & (res$n_secrete == 0 | res$n_uptake == 0))
print(paste0("Reporter metabolites (p_corrected < 0.05): ", sum(rm.all$all.met$p_corrected < 0.05, na.rm = T),
             " in the entire network, ", sum(rm.uptake$all.met$p_corrected < 0.05, na.rm = T), " in the uptake subnetwork"))
print(paste0("Metabolites linked to all species in one direction only (adjusted Z = NA, P = 1): ",
             one_direction_all(rm.all$all.met), " in the entire network, ",
             one_direction_all(rm.uptake$all.met), " in the uptake subnetwork"))
print(paste0("Results are stored in ", output_dir))
