#hgt_domain_counter_and_plotter

library(dplyr)
library(ggplot2)

#Before this script, the reubwn HGT output candidate files were cut -f 1,12 > alienspecies.tab
#Load files and count gene number per hgt candidate (Domain;Kingdom;Species)
Pdetritophagus_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/pdetritophagus_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(Pdetritophagus_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

Pdetritophagus_HGT_species_count_table <- Pdetritophagus_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(Pdetritophagus_HGT_species_count_table)

LC92_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/LC92_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(LC92_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

LC92_HGT_species_count_table <- LC92_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(LC92_HGT_species_count_table)

JU765_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/JU765_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(JU765_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

JU765_HGT_species_count_table <- JU765_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(JU765_HGT_species_count_table)

psuperbus_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/psuperbus_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(psuperbus_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

psuperbus_HGT_species_count_table <- psuperbus_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(psuperbus_HGT_species_count_table)

ES5_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/ES5_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(ES5_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

ES5_HGT_species_count_table <- ES5_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(ES5_HGT_species_count_table)

PS1579_HapA_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/PS1579_HapA_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(PS1579_HapA_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

PS1579_HapA_HGT_species_count_table <- PS1579_HapA_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(PS1579_HapA_HGT_species_count_table)

PS1579_HapB_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/PS1579_HapB_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(PS1579_HapB_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

PS1579_HapB_HGT_species_count_table <- PS1579_HapB_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(PS1579_HapB_HGT_species_count_table)

PS1579_HapC_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/PS1579_HapC_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(PS1579_HapC_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

PS1579_HapC_HGT_species_count_table <- PS1579_HapC_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(PS1579_HapC_HGT_species_count_table)

ps1159_nonem_reubwn_alienspecies <- read.delim("~/alienspecies_tables/ps1159_nonem_reubwn_alienspecies.tab", header=FALSE, comment.char="#")
colnames(ps1159_nonem_reubwn_alienspecies) <- c("Gene ID", "HGT Species")

ps1159_HGT_species_count_table <- ps1159_nonem_reubwn_alienspecies %>%
  group_by(HGT_Species = `HGT Species`) %>%
  summarize(Count = n(), 
            Gene_IDs = paste(`Gene ID`, collapse = ";")) %>%
  ungroup() %>%
  select(Count, HGT_Species, Gene_IDs)

View(ps1159_HGT_species_count_table)

#Combining all tables into one

Pdetritophagus_HGT_species_count_table$species <- "P. detritophagus"
LC92_HGT_species_count_table$species <- "LC92"
JU765_HGT_species_count_table$species <- "JU765"
psuperbus_HGT_species_count_table$species <- "P. superbus"
ES5_HGT_species_count_table$species <- "ES5"
PS1579_HapA_HGT_species_count_table$species <- "PS1579_HapA"
PS1579_HapB_HGT_species_count_table$species <- "PS1579_HapB"
PS1579_HapC_HGT_species_count_table$species <- "PS1579_HapC"
ps1159_HGT_species_count_table$species <- "PS1159"

combined_all <- rbind(Pdetritophagus_HGT_species_count_table, LC92_HGT_species_count_table, JU765_HGT_species_count_table, psuperbus_HGT_species_count_table, ES5_HGT_species_count_table, PS1579_HapA_HGT_species_count_table, PS1579_HapB_HGT_species_count_table, PS1579_HapC_HGT_species_count_table, ps1159_HGT_species_count_table)

#counting types of domains or kingdoms
combined_all$Domain <- sapply(strsplit(combined_all$HGT_Species, ";"), `[`, 1)
combined_all$Kingdom <- sapply(strsplit(combined_all$HGT_Species, ";"), `[`, 2)

#stacked histogram

# Summarizing counts by species and domain to get number for bar segments
combined_summary <- combined_all %>%
  group_by(species, Domain) %>%
  summarize(Total_Count = sum(Count), .groups = 'drop')

# species order for the plot
combined_summary$species <- factor(combined_summary$species, levels = c("P. detritophagus" , "LC92" , "JU765", "P. superbus", "ES5", "PS1579_HapA", "PS1579_HapB", "PS1579_HapC", "PS1159"))

# domain order for the plot
combined_summary$Domain <- factor(combined_summary$Domain, levels = c( "Viruses", "Bacteria", "Eukaryota",  "Archaea", "undef"))

# Choosing color
color_mapping <- c("Archaea" = "#cc0000", "Eukaryota" = "#232ec2", "Viruses" = "yellow", "Bacteria" = "#d65510", "undef" = "#7b7b7b")

# Generating histogram
ggplot(combined_summary, aes(x = species, y = Total_Count, fill = Domain)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) + # bar width
  scale_fill_manual(values = color_mapping) +
  geom_text(aes(label = Total_Count), position = position_stack(vjust = 0.5), color = "black", size = 3.5) + # Add count number in each segment
  theme_minimal() +
  labs(title = "Count of HGT Candidates in Panagrolaimidae by Domain",
       x = "Species",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), # Improve label readability
        plot.margin = margin(10, 20, 10, 20), # Adjust plot margins for more space around
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill = "#e8e8e8", colour = "#e8e8e8"), # Change background color
        panel.grid.major = element_line(colour = "white"), # Change grid color
        panel.grid.minor = element_line(colour = "white"))
