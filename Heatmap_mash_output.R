mash_table1 <- read.delim("~/mash_table1.tsv")
mash_table2 <- read.delim("~/mash_table2.tsv")
mash_table3 <- read.delim("~/mash_table3.tsv")

library(ggplot2)
library(reshape2)
library(dplyr)

# Inverting Mash_Distance to show closeness as higher score
mash_table1$Inverted_Score <- 1 - mash_table1$Mash_Distance
mash_table2$Inverted_Score <- 1 - mash_table2$Mash_Distance
mash_table3$Inverted_Score <- 1 - mash_table3$Mash_Distance

data_melted1 <- melt(mash_table1, id.vars = c("Reference", "Query"), value.name = "Inverted_Score")
data_for_heatmap <- filter(data_melted1, variable == "Inverted_Score")
distance1 <- read.delim("~/mash_table1.tsv")

data_melted2 <- melt(mash_table2, id.vars = c("Reference", "Query"), value.name = "Inverted_Score")
data_for_heatmap2 <- filter(data_melted2, variable == "Inverted_Score")
distance2 <- read.delim("~/mash_table2.tsv")

data_melted3 <- melt(mash_table3, id.vars = c("Reference", "Query"), value.name = "Inverted_Score")
data_for_heatmap3 <- filter(data_melted3, variable == "Inverted_Score")
distance3 <- read.delim("~/mash_table3.tsv")

data_complete <- merge(data_for_heatmap, distance1, by = c("Reference", "Query"), all = TRUE)

data_complete2 <- merge(data_for_heatmap2, distance2, by = c("Reference", "Query"), all = TRUE)

data_complete3 <- merge(data_for_heatmap3, distance3, by = c("Reference", "Query"), all = TRUE)

complete_grid <- expand.grid(Reference = unique(data_complete$Reference),
                             Query = unique(data_complete$Query))

complete_grid2 <- expand.grid(Reference = unique(data_complete2$Reference),
                              Query = unique(data_complete2$Query))

complete_grid3 <- expand.grid(Reference = unique(data_complete3$Reference),
                              Query = unique(data_complete3$Query))

data_with_NAs <- merge(complete_grid, data_complete, 
                       by = c("Reference", "Query"), 
                       all = TRUE)

data_with_NAs2 <- merge(complete_grid2, data_complete2,
                        by = c ("Reference", "Query"),
                        all = TRUE)

data_with_NAs3 <- merge(complete_grid3, data_complete3,
                       by = c ("Reference", "Query"),
                       all = TRUE)

ggplot(data_with_NAs, aes(x = Query, y = reorder(Reference, -as.numeric(Reference)), fill = Mash_Distance)) +
  geom_tile(color = "grey", linewidth = 0.2) + 
  geom_text(aes(label = sprintf("%.3f", Mash_Distance)), size = 2.5, color = "black", na.rm = TRUE) +
  scale_fill_gradient(low = "red", high = "pink", na.value = "white", name = "Mash\nDistance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)

ggplot(data_with_NAs2, aes(x = Query, y = reorder(Reference, -as.numeric(Reference)), fill = Mash_Distance)) +
  geom_tile(color = "grey", linewidth = 0.2) + 
  geom_text(aes(label = sprintf("%.3f", Mash_Distance)), size = 2.5, color = "black", na.rm = TRUE) +
  scale_fill_gradient(low = "red", high = "#ffe6ea", na.value = "white", name = "Mash\nDistance", limits = c(0.10, 1.00), breaks = seq(0.1, 1.00, 0.2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)

ggplot(data_with_NAs3, aes(x = Query, y = reorder(Reference, -as.numeric(Reference)), fill = Mash_Distance)) +
  geom_tile(color = "grey", linewidth = 0.2) + 
  geom_text(aes(label = sprintf("%.3f", Mash_Distance)), size = 2.5, color = "black", na.rm = TRUE) +
  scale_fill_gradient(low = "red", high = "pink", na.value = "white", name = "Mash\nDistance" ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)
