library(readxl)
library(dplyr)
library(here)
library(stringr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(stats)
library(cowplot)
library(patchwork)

#Read in primer efficiency data. calculate coverted primer efficiency 
efficiency_data <- readxl::read_excel(here::here("xlsx_folder/efficiency/efficiency.xlsx"))

efficiency_data$cpe <- efficiency_data$efficiency/100 + 1

# Set the file names for each file
file_list <- list.files("xlsx_folder", pattern = "*.xlsx", full.names = TRUE)

# Create an empty list to hold the data frames
data_list <- list()

# Iterate through each file in the list and read in the data
for(file in file_list) {
  data <- read_excel(file, skip = 24)
  data$date <- as.numeric(gsub(".xlsx", "", basename(file), fixed = TRUE))
  data_list[[basename(file)]] <- data
}

# Combine the data frames into a single data frame
all_data <- bind_rows(data_list, .id = "file_name")

#standardize gene names (all lower case) in Target

all_data <- all_data %>% 
  mutate(Target = tolower(Target)) %>%
  mutate(Sample = toupper(Sample))

# Remove numbers and spaces from a column in the data frame
all_data$Sample <- str_replace_all(all_data$Sample, "[0-9 ]+", "")

#calculate DMSO_CT
DMSO_CT <- all_data %>%
  filter(Sample == "DMSO") %>%
  group_by(Target, Sample) %>%
  summarize(DMSO_CT = mean(`Cq Mean`, na.rm = TRUE), .groups = "drop")

#join DMSO_Cq
all_data <- left_join(all_data, DMSO_CT %>% select(Target,DMSO_CT), by = c("Target"))

#calculate difference in Cq from DMSO average
all_data <- all_data %>% mutate(delta_CT = DMSO_CT - all_data$`Cq Mean`)

#add efficiency data to the df
# Select only the "Target" and "cpe" columns from "efficiency_data"
efficiency_subset <- efficiency_data[c("Target", "cpe", "interest")]

# Merge "all_data" and "efficiency_subset" on the "Target" column
all_data <- merge(all_data, efficiency_subset, by = "Target", all.x = TRUE)


#calculate RQ
all_data <- all_data %>% 
  mutate(RQ = cpe^delta_CT)

# Filter data frame by "Cq SD" column and eliminate values greater than 0.2
all_data<- all_data %>%
  filter(all_data$`Cq SD` <= 0.5)

#calculate Ref RQ
# Calculate the mean "RQ" value for each unique combination of "date" and "Sample" columns in the "efficiency_data" data frame
# Create a new data frame with only rows that have "housekeeping" in the "interest" column
housekeeping_data <- all_data %>%
  filter(interest == "housekeeping")

# Calculate the mean "RQ" value for each unique combination of "date" and "Sample" columns
housekeeping_data <- housekeeping_data %>%
  group_by(date, Sample) %>%
  summarize(mean_RQ = mean(RQ, na.rm = TRUE))

# Join the "RQ_Ref" values from "ref_data" to "all_data" based on the "date" and "Sample" columns
all_data <- all_data %>%
  left_join(housekeeping_data, by = c("date", "Sample"))


# Create the new "GE" column by dividing "RQ" by "RQ_Ref"
all_data <- all_data %>% 
  mutate(GE = RQ/mean_RQ)

key<- read_xlsx(here("key.xlsx"))

all_data <- all_data %>%
  mutate(date = gsub("[[:punct:]]|[a-zA-Z]|^0+", "", file_name))

key_renamed <- key %>%
  rename(date = `file name`)

merged_data <- merge(all_data, key_renamed, by = "date")
merged_data$Sample <- gsub("\\s+", "", merged_data$Sample)


key<- read_xlsx(here("key.xlsx"))

combined_data <- all_data %>%
  mutate(date = gsub("[[:punct:]]|[a-zA-Z]|^0+", "", file_name))

key_renamed <- key %>%
  rename(date = `file name`)

merged_data <- merge(combined_data, key_renamed, by = "date")
merged_data$Sample <- gsub("\\s+", "", merged_data$Sample)

# Convert 'Target' and 'Sample' columns to factors
merged_data$Target <- as.factor(merged_data$Target)
merged_data$Sample <- as.factor(merged_data$Sample)

# Remove rows with missing 'GE' values
merged_data <- merged_data[!is.na(merged_data$GE), ]

## Create 'Pretreatment' column
merged_data$Pretreatment <- ifelse(merged_data$Sample %in% c("DMSO", "NOC"), "DMSO", "DLKI")

# Create 'Treatment' column
merged_data$Treatment <- ifelse(merged_data$Sample %in% c("DMSO", "DLKI"), "DMSO", "NOC")

# Convert 'Pretreatment' and 'Treatment' columns to factors
merged_data$Pretreatment <- as.factor(merged_data$Pretreatment)
merged_data$Treatment <- as.factor(merged_data$Treatment)

# Perform two-way ANOVA for each target in merged_data, excluding specific targets
results <- list()

unique_targets <- setdiff(unique(merged_data$Target), c("flnc", "phldb2", "sipa1l1", "tes"))

for (target in unique_targets) {
  df <- merged_data[merged_data$Target == target, ]
  
  # Convert 'Treatment' and 'Pretreatment' columns to factors with 2 levels
  df$Treatment <- factor(df$Treatment, levels = c("DMSO", "NOC"))
  df$Pretreatment <- factor(df$Pretreatment, levels = c("DMSO", "DLKI"))
  
  # Perform two-way ANOVA for GE using Treatment and Pretreatment
  anova_result <- aov(GE ~ Treatment * Pretreatment, data = df)
  
  # Associate ANOVA result with Target
  results[[target]] <- anova_result
}

# Print the ANOVA results for each Target
for (target in names(results)) {
  anova_result <- results[[target]]
  
  cat("ANOVA for Target:", target, "\n")
  print(summary(anova_result))
  cat("\n")
}

# Subset 1: Targets zyx, itga9, and stmn4 (noc sig)
subset1 <- subset(merged_data, Target %in% c("zyx", "itga9", "stmn4"))

# Subset 2: Targets sh2d3c, sh3, palld, flnb, and tpm4 (DLK inhibitor effect)
subset2 <- subset(merged_data, Target %in% c("sh2d3c", "sh3", "flrt3", "gab1", "palld", "flnb", "tpm4"))

# Reorder the levels of Pretreatment
subset1$Pretreatment <- factor(subset1$Pretreatment, levels = c("DMSO", "DLKI"))
subset1$Target <- factor(subset1$Target, levels = c("stmn4", "zyx", "itga9"))

# Plot data subset1
plot <- ggplot(subset1, aes(x = Pretreatment, y = GE, fill = Treatment, color = Treatment)) +
  geom_boxplot() +
  facet_grid(Target ~ ., switch = "y") +
  theme_bw() +
  labs(x = "Pretreatment", y = "GE") +
  theme(plot.margin = margin(1, 1, 30, 1),
        legend.position = "bottom")


# Print the vertical plot
print(plot)

