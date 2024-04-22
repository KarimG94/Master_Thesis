library(ggplot2)
library(dplyr)
library(tidyr)

# Loading the gff3 files
HapA_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/PS1579.HapA.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
HapB_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/PS1579.HapB.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
HapC_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/PS1579.HapC.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
ES5_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/ES5.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
Psuperbus_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/psuperbus.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
Pdetritophagus_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/pdetritophagus.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
LC92_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/LC92.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
JU765_gff3_data <- read.table("/home/karim/corrected_and_named_protein_files/gff_files/correct_files/JU765.braker.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#", col.names = c("SeqID", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))

# Extracting gene ID for grouping
HapA_gff3_data$GeneID <- sapply(strsplit(as.character(HapA_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE) #replace by Parent= for exon length per gene
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

HapB_gff3_data$GeneID <- sapply(strsplit(as.character(HapB_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE)
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

HapC_gff3_data$GeneID <- sapply(strsplit(as.character(HapC_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE)
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

ES5_gff3_data$GeneID <- sapply(strsplit(as.character(ES5_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE)
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

Psuperbus_gff3_data$GeneID <- sapply(strsplit(as.character(Psuperbus_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE)
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

LC92_gff3_data$GeneID <- sapply(strsplit(as.character(LC92_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE)
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

JU765_gff3_data$GeneID <- sapply(strsplit(as.character(JU765_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE)
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

Pdetritophagus_gff3_data$GeneID <- sapply(strsplit(as.character(Pdetritophagus_gff3_data$Attributes), ";"), function(x) {
  parent <- grep("ID=", x, value = TRUE)
  if(length(parent) > 0){
    return(gsub("ID=", "", parent))
  } else {
    return(NA)
  }
})

# Filtering exons and introns
HapA_exons_introns <- HapA_gff3_data %>% filter(Feature %in% c("exon", "intron"))
HapB_exons_introns <- HapB_gff3_data %>% filter(Feature %in% c("exon", "intron"))
HapC_exons_introns <- HapC_gff3_data %>% filter(Feature %in% c("exon", "intron"))
ES5_exons_introns <- ES5_gff3_data %>% filter(Feature %in% c("exon", "intron"))
Psuperbus_exons_introns <- Psuperbus_gff3_data %>% filter(Feature %in% c("exon", "intron"))
Pdetritophagus_exons_introns <- Pdetritophagus_gff3_data %>% filter(Feature %in% c("exon", "intron"))
LC92_exons_introns <- LC92_gff3_data %>% filter(Feature %in% c("exon", "intron"))
JU765_exons_introns <- JU765_gff3_data %>% filter(Feature %in% c("exon", "intron"))

# Calculating lengths
HapA_exons_introns$Length <- with(HapA_exons_introns, End - Start + 1)
HapB_exons_introns$Length <- with(HapB_exons_introns, End - Start + 1)
HapC_exons_introns$Length <- with(HapC_exons_introns, End - Start + 1)
ES5_exons_introns$Length <- with(ES5_exons_introns, End - Start + 1)
Psuperbus_exons_introns$Length <- with(Psuperbus_exons_introns, End - Start + 1)
LC92_exons_introns$Length <- with(LC92_exons_introns, End - Start + 1)
JU765_exons_introns$Length <- with(JU765_exons_introns, End - Start + 1)
Pdetritophagus_exons_introns$Length <- with(Pdetritophagus_exons_introns, End - Start + 1)

# Add up the lengths of exons and introns for each gene
HapA_lengths_by_gene <- HapA_exons_introns %>%
  group_by(GeneID, Feature) %>% #use ParentID for length of exons per gene
  summarise(TotalLength = sum(Length), .groups = 'drop')

HapB_lengths_by_gene <- HapB_exons_introns %>%
  group_by(GeneID, Feature) %>% 
  summarise(TotalLength = sum(Length), .groups = 'drop')

HapC_lengths_by_gene <- HapC_exons_introns %>%
  group_by(GeneID, Feature) %>% 
  summarise(TotalLength = sum(Length), .groups = 'drop')

ES5_lengths_by_gene <- ES5_exons_introns %>%
  group_by(GeneID, Feature) %>% 
  summarise(TotalLength = sum(Length), .groups = 'drop')

Psuperbus_lengths_by_gene <- Psuperbus_exons_introns %>%
  group_by(GeneID, Feature) %>% 
  summarise(TotalLength = sum(Length), .groups = 'drop')

LC92_lengths_by_gene <- LC92_exons_introns %>%
  group_by(GeneID, Feature) %>% 
  summarise(TotalLength = sum(Length), .groups = 'drop')

JU765_lengths_by_gene <- JU765_exons_introns %>%
  group_by(GeneID, Feature) %>% 
  summarise(TotalLength = sum(Length), .groups = 'drop')

Pdetritophagus_lengths_by_gene <- Pdetritophagus_exons_introns %>%
  group_by(GeneID, Feature) %>% 
  summarise(TotalLength = sum(Length), .groups = 'drop')

# Pivoting
HapA_lengths_wide <- HapA_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

HapB_lengths_wide <- HapB_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

HapC_lengths_wide <- HapC_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

ES5_lengths_wide <- ES5_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

Psuperbus_lengths_wide <- Psuperbus_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

LC92_lengths_wide <- LC92_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

JU765_lengths_wide <- JU765_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

Pdetritophagus_lengths_wide <- Pdetritophagus_lengths_by_gene %>%
  pivot_wider(names_from = Feature, values_from = TotalLength, values_fill = list(TotalLength = 0))

# Calculate average lengths on all genes
HapA_average_lengths <- HapA_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

HapB_average_lengths <- HapB_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

HapC_average_lengths <- HapC_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

ES5_average_lengths <- ES5_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

Psuperbus_average_lengths <- Psuperbus_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

LC92_average_lengths <- LC92_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

JU765_average_lengths <- JU765_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

Pdetritophagus_average_lengths <- Pdetritophagus_lengths_wide %>%
  summarise(AverageExonLength = mean(exon, na.rm = TRUE),
            AverageIntronLength = mean(intron, na.rm = TRUE))

# Finding the end position of the last gene on each chromosome

HapA_last_gene_positions <- HapA_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

HapB_last_gene_positions <- HapB_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

HapC_last_gene_positions <- HapC_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

ES5_last_gene_positions <- ES5_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

Psuperbus_last_gene_positions <- Psuperbus_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

LC92_last_gene_positions <- LC92_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

JU765_last_gene_positions <- JU765_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

Pdetritophagus_last_gene_positions <- Pdetritophagus_gff3_data %>%
  filter(Feature == "gene") %>%
  group_by(SeqID) %>%
  summarise(ChromosomeEnd = max(End), .groups = 'drop')

# Getting total genome length by adding all chromosome lengths
HapA_total_genome_length <- sum(HapA_last_gene_positions$ChromosomeEnd)
HapB_total_genome_length <- sum(HapB_last_gene_positions$ChromosomeEnd)
HapC_total_genome_length <- sum(HapC_last_gene_positions$ChromosomeEnd)
ES5_total_genome_length <- sum(ES5_last_gene_positions$ChromosomeEnd)
Psuperbus_total_genome_length <- sum(Psuperbus_last_gene_positions$ChromosomeEnd)
LC92_total_genome_length <- sum(LC92_last_gene_positions$ChromosomeEnd)
JU765_total_genome_length <- sum(JU765_last_gene_positions$ChromosomeEnd)
Pdetritophagus_total_genome_length <- sum(Pdetritophagus_last_gene_positions$ChromosomeEnd)


# Calculating average intergenic space
# First gene rows are filtered

HapA_sorted <- HapA_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

HapB_sorted <- HapB_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

HapC_sorted <- HapC_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

ES5_sorted <- ES5_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

Psuperbus_sorted <- Psuperbus_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

LC92_sorted <- LC92_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

JU765_sorted <- JU765_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

Pdetritophagus_sorted <- Pdetritophagus_gff3_data %>% 
  filter(Feature == "gene") %>% 
  arrange(SeqID, Start)

#Calculate space
HapA_sorted <- HapA_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start), # Get of next gene on chromosome
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

HapB_sorted <- HapB_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start),
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

HapC_sorted <- HapC_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start),
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

ES5_sorted <- ES5_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start),
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

Psuperbus_sorted <- Psuperbus_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start),
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

LC92_sorted <- LC92_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start),
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

JU765_sorted <- JU765_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start),
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

Pdetritophagus_sorted <- Pdetritophagus_sorted %>%
  group_by(SeqID) %>%
  mutate(NextStart = lead(Start),
         IntergenicSpace = ifelse(row_number() == n(), NA, NextStart - End - 1)) %>%
  ungroup()

# Calculating the average intergenic space

HapA_avg_intergenic <- mean(HapA_sorted$IntergenicSpace, na.rm = TRUE)
HapB_avg_intergenic <- mean(HapB_sorted$IntergenicSpace, na.rm = TRUE)
HapC_avg_intergenic <- mean(HapC_sorted$IntergenicSpace, na.rm = TRUE)
ES5_avg_intergenic <- mean(ES5_sorted$IntergenicSpace, na.rm = TRUE)
Psuperbus_avg_intergenic <- mean(Psuperbus_sorted$IntergenicSpace, na.rm = TRUE)
LC92_avg_intergenic <- mean(LC92_sorted$IntergenicSpace, na.rm = TRUE)
JU765_avg_intergenic <- mean(JU765_sorted$IntergenicSpace, na.rm = TRUE)
Pdetritophagus_avg_intergenic <- mean(Pdetritophagus_sorted$IntergenicSpace, na.rm = TRUE)

# Preparing data frames for plotting
plot_data <- data.frame(
  HapA_FeatureType = c("Exon", "Intron"),
  HapA_AverageLength = c(HapA_average_lengths$AverageExonLength, HapA_average_lengths$AverageIntronLength),
  HapA_GenomeLength = HapA_total_genome_length,
  HapB_FeatureType = c("Exon", "Intron"),
  HapB_AverageLength = c(HapB_average_lengths$AverageExonLength, HapB_average_lengths$AverageIntronLength),
  HapB_GenomeLength = HapB_total_genome_length,
  HapC_FeatureType = c("Exon", "Intron"),
  HapC_AverageLength = c(HapC_average_lengths$AverageExonLength, HapC_average_lengths$AverageIntronLength),
  HapC_GenomeLength = HapC_total_genome_length,
  ES5_FeatureType = c("Exon", "Intron"),
  ES5_AverageLength = c(ES5_average_lengths$AverageExonLength, ES5_average_lengths$AverageIntronLength),
  ES5_GenomeLength = ES5_total_genome_length,
  Psuperbus_FeatureType = c("Exon", "Intron"),
  Psuperbus_AverageLength = c(Psuperbus_average_lengths$AverageExonLength, Psuperbus_average_lengths$AverageIntronLength),
  Psuperbus_GenomeLength = Psuperbus_total_genome_length,
  LC92_FeatureType = c("Exon", "Intron"),
  LC92_AverageLength = c(LC92_average_lengths$AverageExonLength, LC92_average_lengths$AverageIntronLength),
  LC92_GenomeLength = LC92_total_genome_length,
  JU765_FeatureType = c("Exon", "Intron"),
  JU765_AverageLength = c(JU765_average_lengths$AverageExonLength, JU765_average_lengths$AverageIntronLength),
  JU765_GenomeLength = JU765_total_genome_length,
  Pdetritophagus_FeatureType = c("Exon", "Intron"),
  Pdetritophagus_AverageLength = c(Pdetritophagus_average_lengths$AverageExonLength, Pdetritophagus_average_lengths$AverageIntronLength),
  Pdetritophagus_GenomeLength = Pdetritophagus_total_genome_length
)

corrected_data <- plot_data %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Species", ".value"),
    names_sep = "_",
  )

corrected_data$GenomeLength <- (corrected_data$GenomeLength)/10^6
corrected_data$Species <- c("PS1579 HapA" , "PS1579 HapB" , "PS1579 HapC", "ES5" , "P. superbus" , "LC92", "JU765", "P. detritophagus" , "PS1579 HapA" , "PS1579 HapB" , "PS1579 HapC", "ES5" , "P. superbus" , "LC92", "JU765", "P. detritophagus")

# Filtering the data for Exons and Introns
exon_data <- corrected_data %>% filter(FeatureType == "Exon")
intron_data <- corrected_data %>% filter(FeatureType == "Intron")

# Intergenic data frame

intergenic_plot_data <- data.frame(
  Species = c("PS1579 HapA", "PS1579 HapB", "PS1579 HapC", "ES5", "P. superbus", "LC92", "JU765", "P. detritophagus"),
  AverageIntergenicSpace = c(HapA_avg_intergenic, HapB_avg_intergenic, HapC_avg_intergenic, ES5_avg_intergenic, Psuperbus_avg_intergenic, LC92_avg_intergenic, JU765_avg_intergenic, Pdetritophagus_avg_intergenic),
  GenomeLength = c(HapA_total_genome_length, HapB_total_genome_length, HapC_total_genome_length, ES5_total_genome_length, Psuperbus_total_genome_length, LC92_total_genome_length, JU765_total_genome_length, Pdetritophagus_total_genome_length) / 10^6
)

# Calculating Spearman correlation coeffcients between exon/intron/intergenic lengths and genome length

library(Hmisc)

exon_spearman <- cor.test(exon_data$AverageLength, exon_data$GenomeLength, method = "spearman")
print(exon_spearman)
exon_spearman_rho <- -0.76
exon_spearman_pvalue <- 0.037


intron_spearman <- cor.test(intron_data$AverageLength, intron_data$GenomeLength, method = "spearman")
print(intron_spearman)
intron_spearman_rho <- 0.26
intron_spearman_pvalue <- 0.54

intergenic_spearman <- cor.test(intergenic_plot_data$AverageIntergenicSpace, intergenic_plot_data$GenomeLength, method = "spearman")
intergenic_spearman_rho <- 0.79
intergenic_spearman_pvalue <- 0.028

# Plotting
# Colors for species
species_colors <- c("PS1579 HapA" = "#0A75AD", "PS1579 HapB" = "#e51607", "PS1579 HapC" = "#668F25", 
                    "ES5" = "#ceda09", "P. superbus" = "#080870", "LC92" = "orange", 
                    "JU765" = "brown", "P. detritophagus" = "#FF7373")

# Plot for Exons

ggplot(exon_data, aes(x = GenomeLength, y = AverageLength)) +
  geom_point(aes(fill = Species), size = 5, shape = 21, color = "black") +
  scale_fill_manual(values = species_colors) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Total Genome Length (Mbp)", y = "Average Exon Length (bp)", 
       title = "Exon Length to Genome Length", fill = "Species") +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(color = "black"))) + #black border
  theme(legend.text = element_text(face = "italic"))+
  annotate("text", x = 108, y = 175, label = paste("ρ =", format(exon_spearman_rho, digits = 3), 
                                                   "\np-value =", format(exon_spearman_pvalue, digits = 4)), 
           size = 3.5, color = "black")



# Plot for Introns

ggplot(intron_data, aes(x = GenomeLength, y = AverageLength)) +
  geom_point(aes(fill = Species), size = 5, shape = 21, color = "black") +
  scale_fill_manual(values = species_colors) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Total Genome Length (Mbp)", y = "Average Intron Length (bp)", 
       title = "Intron Length to Genome Length", fill = "Species") +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(color = "black"))) +
  theme(legend.text = element_text(face = "italic")) +
  annotate("text", x = 108, y = 80, label = paste("ρ =", format(intron_spearman_rho, digits = 3), 
                                                   "\np-value =", format(intron_spearman_pvalue, digits = 3)), 
           size = 3.5, color = "black")

#Plot for intergenic

ggplot(intergenic_plot_data, aes(x = GenomeLength, y = AverageIntergenicSpace)) +
  geom_point(aes(fill = Species), size = 5, shape = 21, color = "black") +
  scale_fill_manual(values = species_colors) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Total Genome Length (Mbp)", y = "Average Intergenic Space (bp)", 
       title = "Intergenic Space to Genome Length", fill = "Species") +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(color = "black"))) +
  theme(legend.text = element_text(face = "italic")) +
  annotate("text", x = 108, y = 1500, label = paste("ρ =", format(intergenic_spearman_rho, digits = 3), 
                                                  "\np-value =", format(intergenic_spearman_pvalue, digits = 4)), 
           size = 3.5, color = "black")

# Calculating total length of all exons, introns and intergenic and saving as table

#HapA
HapA_total_exons_length <- (sum(HapA_exons_introns$Length[HapA_exons_introns$Feature == "exon"]))/10^6
HapA_total_introns_length <- (sum(HapA_exons_introns$Length[HapA_exons_introns$Feature == "intron"]))/10^6
HapA_total_intergenic_length <- (sum(HapA_sorted$IntergenicSpace, na.rm = TRUE))/10^6

#HapB
HapB_total_exons_length <- (sum(HapB_exons_introns$Length[HapB_exons_introns$Feature == "exon"]))/10^6
HapB_total_introns_length <- (sum(HapB_exons_introns$Length[HapB_exons_introns$Feature == "intron"]))/10^6
HapB_total_intergenic_length <- (sum(HapB_sorted$IntergenicSpace, na.rm = TRUE))/10^6

#Hapc
HapC_total_exons_length <- (sum(HapC_exons_introns$Length[HapC_exons_introns$Feature == "exon"]))/10^6
HapC_total_introns_length <- (sum(HapC_exons_introns$Length[HapC_exons_introns$Feature == "intron"]))/10^6
HapC_total_intergenic_length <- (sum(HapC_sorted$IntergenicSpace, na.rm = TRUE))/10^6

#ES5
ES5_total_exons_length <- (sum(ES5_exons_introns$Length[ES5_exons_introns$Feature == "exon"]))/10^6
ES5_total_introns_length <- (sum(ES5_exons_introns$Length[ES5_exons_introns$Feature == "intron"]))/10^6
ES5_total_intergenic_length <- (sum(ES5_sorted$IntergenicSpace, na.rm = TRUE))/10^6

#HapA
Psuperbus_total_exons_length <- (sum(Psuperbus_exons_introns$Length[Psuperbus_exons_introns$Feature == "exon"]))/10^6
Psuperbus_total_introns_length <- (sum(Psuperbus_exons_introns$Length[Psuperbus_exons_introns$Feature == "intron"]))/10^6
Psuperbus_total_intergenic_length <- (sum(Psuperbus_sorted$IntergenicSpace, na.rm = TRUE))/10^6

#LC92
LC92_total_exons_length <- (sum(LC92_exons_introns$Length[LC92_exons_introns$Feature == "exon"]))/10^6
LC92_total_introns_length <- (sum(LC92_exons_introns$Length[LC92_exons_introns$Feature == "intron"]))/10^6
LC92_total_intergenic_length <- (sum(LC92_sorted$IntergenicSpace, na.rm = TRUE))/10^6

#JU765
JU765_total_exons_length <- (sum(JU765_exons_introns$Length[JU765_exons_introns$Feature == "exon"]))/10^6
JU765_total_introns_length <- (sum(JU765_exons_introns$Length[JU765_exons_introns$Feature == "intron"]))/10^6
JU765_total_intergenic_length <- (sum(JU765_sorted$IntergenicSpace, na.rm = TRUE))/10^6

#Pdetritophagus
Pdetritophagus_total_exons_length <- (sum(Pdetritophagus_exons_introns$Length[Pdetritophagus_exons_introns$Feature == "exon"]))/10^6
Pdetritophagus_total_introns_length <- (sum(Pdetritophagus_exons_introns$Length[Pdetritophagus_exons_introns$Feature == "intron"]))/10^6
Pdetritophagus_total_intergenic_length <- (sum(Pdetritophagus_sorted$IntergenicSpace, na.rm = TRUE))/10^6

# Combine into a data frame
length_table <- data.frame(
  Species = c("PS1579 HapA", "PS1579 HapB", "PS1579 HapC", "ES5" , "P. superbus", "JU765", "LC92", "P. detritophagus"),
  ExonLength = c(HapA_total_exons_length, HapB_total_exons_length, HapC_total_exons_length, ES5_total_exons_length, Psuperbus_total_exons_length, JU765_total_exons_length, LC92_total_exons_length, Pdetritophagus_total_exons_length),
  IntronLength = c(HapA_total_introns_length, HapB_total_introns_length, HapC_total_introns_length, ES5_total_introns_length, Psuperbus_total_introns_length, JU765_total_introns_length, LC92_total_introns_length, Pdetritophagus_total_introns_length),
  IntergenicSpace = c(HapA_total_intergenic_length, HapB_total_intergenic_length, HapC_total_intergenic_length, ES5_total_intergenic_length, Psuperbus_total_intergenic_length, JU765_total_intergenic_length, LC92_total_intergenic_length, Pdetritophagus_total_intergenic_length)
)

write.table(length_table, file="species_lengths.tsv", sep="\t", row.names=FALSE, quote=FALSE)
