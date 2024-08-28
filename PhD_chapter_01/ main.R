library(rtracklayer)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(FSA)
library(tidyverse)
library(dunn.test)
library(readr)

# Function to process each GFF file
process_gff_file <- function(gff_file, species_name) {
  # Try to read the GFF file and handle any errors
  tryCatch({
    # Read the GFF file
    gff_data <- import(gff_file)
    
    # Convert the GFF data to a data frame
    gff_df <- as.data.frame(gff_data)
    
    # Remove redundant rows from gff_df
    gff_df <- gff_df %>%
      distinct(seqnames, start, end, type, strand, gene, Name, .keep_all = TRUE)
    
    # Get chromosome sizes
    chromosome_sizes <- gff_df %>%
      group_by(seqnames) %>%
      summarize(chromosome_size = max(end)) %>%
      ungroup()
    
    # Initialize an empty data frame to store the results
    results <- data.frame(
      species = character(),
      chromosome = character(),
      chromosome_size = integer(),
      lncRNA = character(),
      gene_ID = character(),
      gene_bank = character(),
      lncRNA_strand = character(),
      lncRNA_start = integer(),
      lncRNA_end = integer(),
      previous_mRNA = character(),
      p_mRNA_strand = character(),
      previous_mRNA_start = integer(),
      previous_mRNA_end = integer(),
      next_mRNA = character(),
      n_mRNA_strand = character(),
      next_mRNA_start = integer(),
      next_mRNA_end = integer(),
      strand_combination = character(),
      lncRNA_p_mRNA = character(),
      lnRNA_n_mRNA = character(),
      lncRNA_p_mRNA_distance = integer(),
      lncRNA_n_mRNA_distance = integer(),
      closest_distance = integer(),
      random = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Find indices of lnc_RNA rows
    lnc_indices <- which(gff_df$type == "lnc_RNA")
    
    for (lnc_idx in lnc_indices) {
      lncRNA_start <- gff_df$start[lnc_idx]
      lncRNA_end <- gff_df$end[lnc_idx]
      
      # Find the previous mRNA row
      previous_mRNA_idx <- tail(which(gff_df$type == "mRNA" & 1:nrow(gff_df) < lnc_idx), 1)
      previous_mRNA_start <- ifelse(length(previous_mRNA_idx) == 0, NA, gff_df$start[previous_mRNA_idx])
      previous_mRNA_end <- ifelse(length(previous_mRNA_idx) == 0, NA, gff_df$end[previous_mRNA_idx])
      p_mRNA_strand_value <- ifelse(length(previous_mRNA_idx) == 0, NA, as.character(gff_df$strand[previous_mRNA_idx]))
      
      # Find the next mRNA row
      next_mRNA_idx <- head(which(gff_df$type == "mRNA" & 1:nrow(gff_df) > lnc_idx), 1)
      next_mRNA_start <- ifelse(length(next_mRNA_idx) == 0, NA, gff_df$start[next_mRNA_idx])
      next_mRNA_end <- ifelse(length(next_mRNA_idx) == 0, NA, gff_df$end[next_mRNA_idx])
      n_mRNA_strand_value <- ifelse(length(next_mRNA_idx) == 0, NA, as.character(gff_df$strand[next_mRNA_idx]))
      
      # Get the lncRNA gene type value
      lncRNA_value <- gff_df$gene[lnc_idx]
      lncRNA_strand_value <- as.character(gff_df$strand[lnc_idx])
      
      # Get chromosome value
      chromosome_value <- gff_df$seqnames[lnc_idx]
      
      # Get chromosome size value
      chromosome_size_value <- chromosome_sizes$chromosome_size[chromosome_sizes$seqnames == chromosome_value]
      
      # Get the mRNA gene type value for previous mRNA
      previous_mRNA_value <- ifelse(length(previous_mRNA_idx) == 0, NA, gff_df$gene[previous_mRNA_idx])
      
      # Get the mRNA gene type value for next mRNA
      next_mRNA_value <- ifelse(length(next_mRNA_idx) == 0, NA, gff_df$gene[next_mRNA_idx])
      
      # Determine the strand_combination value
      strand_combination_value <- paste0(
        ifelse(is.na(lncRNA_strand_value), "", lncRNA_strand_value),
        ifelse(is.na(p_mRNA_strand_value), "", p_mRNA_strand_value),
        ifelse(is.na(n_mRNA_strand_value), "", n_mRNA_strand_value)
      )
      
      # Determine strand combination for lncRNA and previous mRNA
      lncRNA_p_mRNA_value <- paste0(ifelse(is.na(lncRNA_strand_value), "", lncRNA_strand_value), ifelse(is.na(p_mRNA_strand_value), "", p_mRNA_strand_value))
      
      # Determine strand combination for lncRNA and next mRNA
      lncRNA_n_mRNA_value <- paste0(ifelse(is.na(lncRNA_strand_value), "", lncRNA_strand_value), ifelse(is.na(n_mRNA_strand_value), "", n_mRNA_strand_value))
      
      # Determine the distance between lncRNA and previous mRNA
      lncRNA_p_mRNA_distance_value <- ifelse(!is.na(previous_mRNA_end) & lncRNA_start <= previous_mRNA_end, 0, ifelse(!is.na(previous_mRNA_end) & !is.na(lncRNA_start), lncRNA_start - previous_mRNA_end, NA))
      
      # Determine the distance between lncRNA and next mRNA
      lncRNA_n_mRNA_distance_value <- ifelse(!is.na(next_mRNA_start) & lncRNA_end >= next_mRNA_start, 0, ifelse(!is.na(next_mRNA_start) & !is.na(lncRNA_end), lncRNA_end - next_mRNA_start, NA))
      
      # Calculate closest distance with original sign
      closest_distance_value <- ifelse(abs(lncRNA_p_mRNA_distance_value) < abs(lncRNA_n_mRNA_distance_value), lncRNA_p_mRNA_distance_value, lncRNA_n_mRNA_distance_value)
      
      # Add the results to the data frame
      results <- rbind(results, data.frame(
        species = species_name,
        chromosome = chromosome_value,
        chromosome_size = chromosome_size_value,
        lncRNA = lncRNA_value,
        gene_ID = str_remove(lncRNA_value, "LOC"),
        gene_bank = gff_df$Name[lnc_idx],
        lncRNA_strand = lncRNA_strand_value,
        lncRNA_start = lncRNA_start,
        lncRNA_end = lncRNA_end,
        previous_mRNA = previous_mRNA_value,
        p_mRNA_strand = p_mRNA_strand_value,
        previous_mRNA_start = previous_mRNA_start,
        previous_mRNA_end = previous_mRNA_end,
        next_mRNA = next_mRNA_value,
        n_mRNA_strand = n_mRNA_strand_value,
        next_mRNA_start = next_mRNA_start,
        next_mRNA_end = next_mRNA_end,
        strand_combination = strand_combination_value,
        lncRNA_p_mRNA = lncRNA_p_mRNA_value,
        lnRNA_n_mRNA = lncRNA_n_mRNA_value,
        lncRNA_p_mRNA_distance = lncRNA_p_mRNA_distance_value,
        lncRNA_n_mRNA_distance = lncRNA_n_mRNA_distance_value,
        closest_distance = closest_distance_value,
        random = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    # Remove redundant rows based on previous mRNA, keeping the first occurrence
    results <- results[!duplicated(results$previous_mRNA), ]
    
    # Remove NAs from the lncRNA_p_mRNA_distance column and leave blanks
    results$lncRNA_p_mRNA_distance <- ifelse(is.na(results$lncRNA_p_mRNA_distance), "", results$lncRNA_p_mRNA_distance)
    
    # Remove NAs from the lncRNA_n_mRNA_distance column and leave blanks
    results$lncRNA_n_mRNA_distance <- ifelse(is.na(results$lncRNA_n_mRNA_distance), "", results$lncRNA_n_mRNA_distance)
    
    # Remove rows that have NA values in "closest_distance"
    results <- results %>% drop_na(closest_distance)
    
    # Generate random numbers for the whole genome and assign them to the new column "random_whole_genome"
    results$random <- runif(nrow(results), min = -250000, max = 250000)
    
    return(results)
  }, error = function(e) {
    message(paste("Error processing file:", gff_file))
    message(e)
    return(data.frame())
  })
}

# List of GFF files and species names
gff_files <- list(
  "Peromyscus_eremicus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Peromyscus\ eremicus/genomic.gff",
  "Peromyscus_maniculatus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Peromyscus\ maniculatus/GCF_003704035.1_HU_Pman_2.1.3_genomic.gff",
  "Rattus_norvegicus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Rattus\ norvegicus/GCF_036323735.1_GRCr8_genomic.gff",
  "Mus_musculus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Mus\ musculus/GCF_000001635.27_GRCm39_genomic.gff",
  "Perognathus_pacificus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Perognathus\ longimembris\ pacificus/GCF_023159225.1_ASM2315922v1_genomic.gff",
  "Apodemus_sylvaticus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Apodemus\ sylvaticus/GCF_947179515.1_mApoSyl1.1_genomic.gff",
  "Arvicola_amphibius" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Arvicola\ amphibius/GCF_903992535.2_mArvAmp1.2_genomic.gff",
  "Jaculus_jaculus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Jaculus\ jaculus\ /GCF_020740685.1_mJacJac1.mat.Y.cur_genomic.gff",
  "Ochotona_princeps" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Ochotona\ princeps/GCF_030435755.1_mOchPri1.hap1_genomic.gff",
  "Sciurus_carolinensis" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Sciurus\ carolinensis/GCF_902686445.1_mSciCar1.2_genomic.gff",
  "Desmodus_rotundus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Desmodus\ rotundus/GCF_022682495.1_HLdesRot8A_genomic.gff",
  "Homo_sapiens" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Homo\ sapiens/GCF_000001405.40_GRCh38.p14_genomic.gff",
  "Pan_troglodytes" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Pan\ troglodytes/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gff",
  "Pongo_pygmaeus" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Pongo\ pygmaeus/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gff",
  "Macaca_mulatta" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Macaca\ mulatta/GCF_003339765.1_Mmul_10_genomic.gff",
  "Columba_livia" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Columba\ livia/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.gff",
  "Parus_major" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Parus\ major/GCF_001522545.3_Parus_major1.1_genomic.gff",
  "Apis_mellifera" = "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/gff\ files/Apis\ mellifera/GCF_003254395.2_Amel_HAv3.1_genomic.gff"
  
  # Add more species and their corresponding GFF file paths
)

# Initialize an empty data frame to store combined results
combined_results <- data.frame()

# Process each GFF file and combine results
for (species in names(gff_files)) {
  gff_file <- gff_files[[species]]
  species_results <- process_gff_file(gff_file, species)
  combined_results <- rbind(combined_results, species_results)
}

# Keep only the rows that start with "NC"
combined_results <- combined_results %>% filter(str_starts(chromosome, "NC"))

# Save the combined_results data frame as an Excel file
combined <- "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/results/combined_results.xlsx"
write.xlsx(combined_results, combined)

# Print combined results
print(combined_results)

###----------------------------------------------------------------------------------------------------------------###

###----------------------------------------------------------------------------------------------------------------###

### Remove sex chromosomes from the dataset ###
# List of chromosome values to remove
chromosomes_to_remove <- c("NC_067495.1", "NC_052065.1", "NC_088640.1", "NC_088642.1",
                           "NC_071400.1", "NC_000023.11", "NC_000024.10", "NC_059125.1",
                           "NC_027914.1", "NC_041774.1", "NC_000086.8", "NC_000087.8",
                           "NC_080865.1", "NC_080866.1", "NC_072421.2", "NC_072422.2",
                           "NC_031798.1", "NC_031799.1", "NC_081439.1", "NC_056031.1",
                           "NC_072396.2", "NC_072397.2", "NC_086039.1", "NC_086040.1",
                           "NC_062232.1")

# Remove rows with the specified chromosome values
combined_results_filtered <- combined_results %>%
  filter(!chromosome %in% chromosomes_to_remove)

###----------------------------------------------------------------------------------------------------------------###

###----------------------------------------------------------------------------------------------------------------###

### 1) WHOLE GENOME analysis ###
# for either "combined_results_filtered" or "sampled_combined_results"


# Perform Kruskal-Wallis (KW) test (compare whole chromosomes of different species)
kruskal_test_result <- kruskal.test(closest_distance ~ species, data = combined_results_filtered)

# Print KW test result
print(kruskal_test_result)

# Perform post-hoc pairwise comparisons using the sampled data
posthoc_result <- pairwise.wilcox.test(combined_results_filtered$closest_distance, combined_results_filtered$species, p.adjust.method = "bonferroni")

# Print post-hoc test result
print(posthoc_result)


# save pot-hoc test result

# Define the file path
file_path <- "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/results/whole_genome/wilcox_filtered_results_updated.txt"

# Open the file for writing
file_conn <- file(file_path)

# Write the pairwise Wilcoxon test results to the file
writeLines(capture.output(print(posthoc_result)), file_conn)

# Close the file connection
close(file_conn)


# Extract significant results from post-hoc test
significant_results <- as.data.frame(posthoc_result$p.value) %>%
  rownames_to_column("species1") %>%
  pivot_longer(cols = -species1, names_to = "species2", values_to = "p.value") %>%
  filter(p.value < 0.05, !is.na(p.value))

# Save significant results to a CSV file
significant_results_path <- "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/results/whole_genome/significant_wilcox_filtered_results.csv"
write.csv(significant_results, file = significant_results_path, row.names = FALSE)

# Print significant results
print(significant_results)



### Histograms with a normal distribution overlay ###
### Kolmogorov-Smirnov (KS) ###

# Plot histogram for the WHOLE GENOME of each species and perform Kolmogorov-Smirnov (KS) as well

ks_results_file <- "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/results/whole_genome/ks_test_results.txt"
sink(ks_results_file)

species_list <- unique(combined_results_filtered$species)

for (species in species_list) {
  species_data <- combined_results_filtered %>% filter(species == !!species)
  
  # Perform Kolmogorov-Smirnov (KS) test (test of normality: whole chromosome vs normal distribution)
  ks_test <- ks.test(species_data$closest_distance, "pnorm", mean = mean(species_data$random, na.rm = TRUE), sd = sd(species_data$random, na.rm = TRUE))
  
  # Print KS test results
  cat("KS test for species:", species, "\n")
  print(ks_test)
  
  p <- ggplot(species_data, aes(x = closest_distance)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "blue", alpha = 0.5) +
    stat_function(fun = dnorm, args = list(mean = mean(species_data$random, na.rm = TRUE), 
                                           sd = sd(species_data$random, na.rm = TRUE)), 
                  color = "red", size = 1) +
    xlim(-250000, 250000) +
    labs(title = paste("Histogram for whole genome of", species),
         x = "Closest Distance",
         y = "Density") +
    theme_minimal()
  
  print(p)
}

sink() # Close the sink


###----------------------------------------------------------------------------------------------------------------###

###----------------------------------------------------------------------------------------------------------------###

### 2) INDIVIDUAL CHROMOSOME analysis ###

### Histograms with a normal distribution overlay ###

# Alternative 1: plot histograms for all species at once
for (species in species_list) {
  species_data <- combined_results_filtered %>% filter(species == !!species)
  chromosome_list <- unique(species_data$chromosome)
  
  for (chromosome in chromosome_list) {
    chromosome_data <- species_data %>% filter(chromosome == !!chromosome)
    
    # Filter out non-finite values
    chromosome_data <- chromosome_data %>% filter(is.finite(closest_distance) & is.finite(random))
    
    if (nrow(chromosome_data) > 1) { # Adjusted condition
      p <- ggplot(chromosome_data, aes(x = closest_distance)) +
        geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "blue", alpha = 0.5) +
        stat_function(fun = dnorm, args = list(mean = mean(chromosome_data$random, na.rm = TRUE), 
                                               sd = sd(chromosome_data$random, na.rm = TRUE)), 
                      color = "red", size = 1) +
        xlim(-250000, 250000) +
        labs(title = paste("Histogram for", chromosome, "of", species),
             x = "Closest Distance",
             y = "Density") +
        theme_minimal()
      
      plot_filename <- paste0("/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/results/histogram_", species, "_", chromosome, ".png")
      ggsave(plot_filename, plot = p)
    } else {
      cat("Not enough data for plotting histogram for chromosome:", chromosome, "of species:", species, "\n")
    }
  }
}


# Alternative 2: plot histograms species by species
species_to_plot <- "Apis_mellifera"    # plug in the name of species

species_data <- combined_results_filtered %>% filter(species == species_to_plot)
chromosome_list <- unique(species_data$chromosome)

for (chromosome in chromosome_list) {
  chromosome_data <- species_data %>% filter(chromosome == !!chromosome)
  
  # Filter out non-finite values
  chromosome_data <- chromosome_data %>% filter(is.finite(closest_distance) & is.finite(random))
  
  p <- ggplot(chromosome_data, aes(x = closest_distance)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "blue", alpha = 0.5) +
    stat_function(fun = dnorm, args = list(mean = mean(chromosome_data$random, na.rm = TRUE), 
                                           sd = sd(chromosome_data$random, na.rm = TRUE)), 
                  color = "red", size = 1) +
    xlim(-250000, 250000) +
    labs(title = paste("Histogram for", chromosome, "of", species_to_plot),
         x = "Closest Distance",
         y = "Density") +
    theme_minimal()
  
  print(p)
}


### Kolmogorov-Smirnov (KS) ###

# Initialize lists to store ks test results
ks_test_results <- list()

# Perform Kolmogorov-Smirnov (KS) analysis
species_list <- unique(combined_results_filtered$species)

for (species in species_list) {
  species_data <- combined_results_filtered %>% filter(species == !!species)
  chromosome_list <- unique(species_data$chromosome)
  
  for (chromosome in chromosome_list) {
    chromosome_data <- species_data %>% filter(chromosome == !!chromosome)
    
    # Filter out non-finite values
    chromosome_data <- chromosome_data %>% filter(is.finite(closest_distance) & is.finite(random))
    
    if (nrow(chromosome_data) > 10) { # Adjusted condition
      # Perform Kolmogorov-Smirnov (KS) test
      ks_test <- ks.test(chromosome_data$closest_distance, chromosome_data$random)
      
      # Store KS test results
      ks_test_results[[paste(species, chromosome, sep = "_")]] <- ks_test
      
      # Print KS test results
      cat("KS test for chromosome:", chromosome, "of species:", species, "\n")
      print(ks_test)
    } else {
      cat("Not enough data for KS test for chromosome:", chromosome, "of species:", species, "\n")
    }
  }
}

# Prepare KS test results for writing to Excel
ks_results_df <- data.frame(
  species_chromosome = names(ks_test_results),
  statistic = sapply(ks_test_results, function(x) x$statistic),
  p_value = sapply(ks_test_results, function(x) x$p.value),
  stringsAsFactors = FALSE
)

# Write KS test results to Excel
write.xlsx(ks_results_df, "ks_chromosome.xlsx")



### Kruskal-Wallis (KW) and Dunn's tests ###

# Initialize lists to store kw and dunn test results
kw_test_results <- list()
dunn_test_results <- list()

# Perform Kruskal-Wallis (KW) and Dunn's tests
for (species in species_list) {
  species_data <- combined_results_filtered %>% filter(species == !!species)
  
  # Perform Kruskal-Wallis (KW) test
  kw_test <- kruskal.test(closest_distance ~ chromosome, data = species_data)
  
  # Store KW test results
  kw_test_results[[species]] <- kw_test
  
  # Print KW test results
  cat("Kruskal-Wallis test for species:", species, "\n")
  print(kw_test)
  
  # Perform Dunn's test for species with significant Kruskal-Wallis results
  if (kw_test$p.value < 0.05) {
    dunn_test <- dunn.test(species_data$closest_distance, species_data$chromosome, method = "bonferroni")
    dunn_test_results[[species]] <- dunn_test
  }
}

# Prepare Kruskal-Wallis test results for writing to Excel
kw_results_df <- data.frame(
  species = names(kw_test_results),
  statistic = sapply(kw_test_results, function(x) x$statistic),
  p_value = sapply(kw_test_results, function(x) x$p.value),
  stringsAsFactors = FALSE
)

# Write KW test results to Excel
write.xlsx(kw_results_df, "kw_chromosome.xlsx")

# Prepare Dunn's test results for writing to Excel
dunn_results_df <- data.frame()

for (species in names(dunn_test_results)) {
  result <- dunn_test_results[[species]]
  species_df <- data.frame(
    species = species,
    comparison = result$comparisons,
    z = result$Z,
    p_value = result$P.adjusted,
    stringsAsFactors = FALSE
  )
  dunn_results_df <- rbind(dunn_results_df, species_df)
}

# Write Dunn's test results to Excel
write.xlsx(dunn_results_df, "dunn_chromosome.xlsx")

###----------------------------------------------------------------------------------------------------------------###

###----------------------------------------------------------------------------------------------------------------###

### Distribution for single chromosomes ###

# these chromosomes belong to Homo sapiens and Mus musculus. They are significantly different according to Dunn's test.
# this block of code tries to find why these few chromosomes behave differently than all others.


# Create the distribution_chromosomes dataframe
distribution_chromosomes <- data.frame(
  chromosome = c(
    "NC_000073.7", "NC_000074.7",
    "NC_000077.7", "NC_000078.7", "NC_000083.7",
    "NC_000001.11", "NC_000002.12"
  ),
  species = c(
    "Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus",
    "Mus musculus",
    "Homo sapiens", "Homo sapiens"
  ),
  stringsAsFactors = FALSE
)

# Filter out the chromosomes for Homo sapiens and Mus musculus that are not in the specified list
additional_chromosomes <- combined_results_filtered %>%
  filter(species %in% c("Homo_sapiens", "Mus_musculus")) %>%
  filter(!chromosome %in% distribution_chromosomes$chromosome) %>%
  select(chromosome, species) %>%
  distinct()

# Add the additional chromosomes to the distribution_chromosomes dataframe
distribution_chromosomes <- bind_rows(distribution_chromosomes, additional_chromosomes)

# Count the number of lncRNAs for each chromosome in combined_results
lncRNA_counts <- combined_results %>%
  group_by(chromosome) %>%
  summarize(lncRNA_count = n()) %>%
  ungroup()

# Add the lncRNA counts to the distribution_chromosomes dataframe
distribution_chromosomes <- left_join(distribution_chromosomes, lncRNA_counts, by = "chromosome")

# Rename the lncRNA_count column to "#lncRNA"
distribution_chromosomes <- distribution_chromosomes %>%
  rename(`#lncRNA` = lncRNA_count)

# Count the number of closest_distance values that are 0 for each chromosome
zero_distance_counts <- combined_results %>%
  filter(closest_distance == 0) %>%
  group_by(chromosome) %>%
  summarize(zero_distance_count = n()) %>%
  ungroup()

# Add the zero_distance counts to the distribution_chromosomes dataframe
distribution_chromosomes <- left_join(distribution_chromosomes, zero_distance_counts, by = "chromosome")

# Calculate the percentage, keep two decimals, and add a percentage sign
distribution_chromosomes <- distribution_chromosomes %>%
  mutate(`dist = 0` = ifelse(is.na(zero_distance_count), NA, 
                             sprintf("%.2f%%", (zero_distance_count / `#lncRNA`) * 100))) %>%
  select(-zero_distance_count)

# Count the number of closest_distance values that are less than 1000 for each chromosome
less_than_1000_counts <- combined_results_filtered %>%
  filter(closest_distance < 1000) %>%
  group_by(chromosome) %>%
  summarize(less_than_1000_count = n()) %>%
  ungroup()

# Add the less_than_1000 counts to the distribution_chromosomes dataframe
distribution_chromosomes <- left_join(distribution_chromosomes, less_than_1000_counts, by = "chromosome")

# Calculate the percentage, keep two decimals, and add a percentage sign
distribution_chromosomes <- distribution_chromosomes %>%
  mutate(`dist < 1000` = ifelse(is.na(less_than_1000_count), NA, 
                                sprintf("%.2f%%", (less_than_1000_count / `#lncRNA`) * 100))) %>%
  select(-less_than_1000_count)

# Function to calculate percentages for a given range
calculate_percentage <- function(data, lower_bound, upper_bound) {
  data %>%
    filter(closest_distance >= lower_bound & closest_distance < upper_bound) %>%
    group_by(chromosome) %>%
    summarize(count = n()) %>%
    ungroup()
}

# Calculate percentages for different ranges and add to the dataframe
ranges <- list(
  "1000-10k" = c(1000, 10000),
  "10k-100k" = c(10000, 100000),
  "100k-250k" = c(100000, 250000),
  "> 250k" = c(250000, Inf)
)

for (range in names(ranges)) {
  lower_bound <- ranges[[range]][1]
  upper_bound <- ranges[[range]][2]
  counts <- calculate_percentage(combined_results_filtered, lower_bound, upper_bound)
  distribution_chromosomes <- left_join(distribution_chromosomes, counts, by = "chromosome")
  distribution_chromosomes <- distribution_chromosomes %>%
    mutate(!!paste0("dist ", range) := ifelse(is.na(count), NA, 
                                              sprintf("%.2f%%", (count / `#lncRNA`) * 100))) %>%
    select(-count)
}

# Ensure the percentage columns are numeric by removing the '%' sign and converting to numeric
percentage_columns <- c("dist = 0", "dist < 1000", "dist 1000-10k", "dist 10k-100k", "dist 100k-250k", "dist > 250k")
distribution_chromosomes[percentage_columns] <- lapply(distribution_chromosomes[percentage_columns], function(x) {
  as.numeric(sub("%", "", x))
})

# Calculate mean and standard deviation for each percentage column
distribution_chromosomes <- distribution_chromosomes %>%
  rowwise() %>%
  mutate(`mean +- SD` = paste0(
    sprintf("%.2f", mean(c_across(all_of(percentage_columns)), na.rm = TRUE)), 
    " ± ", 
    sprintf("%.2f", sd(c_across(all_of(percentage_columns)), na.rm = TRUE))
  ))

# Print the updated distribution_chromosomes dataframe
print(distribution_chromosomes)

# Save the distribution_results data frame as an Excel file
dist_path <- "/Users/sahar/Documents/Ph.D./PhD_project/chapter_01/results/dist_results/dist_chromosomes.xlsx"
write.xlsx(distribution_chromosomes, dist_path)



### Violin plot ###

library(ggbeeswarm)
library(reshape2)

# Ensure the percentage columns are numeric by removing the '%' sign and converting to numeric
percentage_columns_1 <- c("dist = 0", "dist < 1000")
percentage_columns_2 <- c("dist 1000-10k", "dist 10k-100k", "dist 100k-250k", "dist > 250k")

distribution_chromosomes[percentage_columns_1] <- lapply(distribution_chromosomes[percentage_columns_1], function(x) {
  as.numeric(sub("%", "", x))
})

distribution_chromosomes[percentage_columns_2] <- lapply(distribution_chromosomes[percentage_columns_2], function(x) {
  as.numeric(sub("%", "", x))
})

# Melt the data into a long format for ggplot2
long_data_1 <- melt(distribution_chromosomes, id.vars = c("species"), measure.vars = percentage_columns_1)
long_data_2 <- melt(distribution_chromosomes, id.vars = c("species"), measure.vars = percentage_columns_2)

# Plot bean plot for dist = 0 and dist < 1000 using ggplot2 and ggbeeswarm
bean_plot_1 <- ggplot(long_data_1, aes(x = variable, y = value, color = species)) +
  geom_violin(trim = FALSE) +  # Violin plot as an alternative to bean plot
  geom_beeswarm(size = 2, alpha = 0.6) +  # Beeswarm plot for individual points
  labs(title = "Bean Plot for Distance = 0 and Distance < 1000",
       x = "Distance Category",
       y = "Percentage",
       color = "Species") +
  theme_minimal()

# Print bean plot 1
print(bean_plot_1)

# Plot bean plot for dist 1000-10k, dist 10k-100k, dist 100k-250k, dist > 250k using ggplot2 and ggbeeswarm
bean_plot_2 <- ggplot(long_data_2, aes(x = variable, y = value, color = species)) +
  geom_violin(trim = FALSE) +  # Violin plot as an alternative to bean plot
  geom_beeswarm(size = 2, alpha = 0.6) +  # Beeswarm plot for individual points
  labs(title = "Bean Plot for Distance 1000-10k, 10k-100k, 100k-250k, > 250k",
       x = "Distance Category",
       y = "Percentage",
       color = "Species") +
  theme_minimal()

# Print bean plot 2
print(bean_plot_2)

###----------------------------------------------------------------------------------------------------------------###
