#### Sex Differnce ####
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')
# Summarize counts by group
table_sex_vhl <- table(tumors$Sex, tumors$VHL)

# Calculate proportions
prop_sex_vhl <- prop.table(table_sex_vhl, margin = 2)

# Print the tables
print(table_sex_vhl)  # Count table
print(round(prop_sex_vhl * 100, 2))  # Percentage table

chisq_test <- chisq.test(table_sex_vhl)
print(chisq_test)

fisher_test <- fisher.test(table_sex_vhl)
print(fisher_test)

library(ggplot2)
ggplot(tumors, aes(x = VHL, fill = Sex)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", title = "Sex Distribution by VHL Status") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal()

ggplot(tumors, aes(x = VHL, fill = Sex)) +
  geom_bar(position = "stack") +
  labs(y = "Count", title = "Sex Counts by VHL Status") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal()


#### Sex Difference Stratified by Stage ####

library(readxl)
library(ggplot2)
library(gridExtra)

# Load the data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')

# Ensure Stage is a factor and ordered
tumors$Stage <- factor(tumors$Stage, levels = c("I", "II", "III", "IV"))

# Stratify analysis by ordered Stage
unique_stages <- levels(tumors$Stage)  # Get the ordered unique stages

# Initialize lists to store plots
proportion_plots <- list()
count_plots <- list()

for (stage in unique_stages) {
  cat("\n### Stage:", stage, "###\n")
  
  # Subset data for the current stage
  tumors_stage <- tumors[tumors$Stage == stage, ]
  
  # Summarize counts by Sex and VHL within this stage
  table_sex_vhl <- table(tumors_stage$Sex, tumors_stage$VHL)
  
  # Skip stages with empty or invalid data
  if (sum(table_sex_vhl) == 0) {
    cat("No valid data for Stage:", stage, "\n")
    next
  }
  
  # Calculate proportions
  prop_sex_vhl <- prop.table(table_sex_vhl, margin = 2)
  
  # Print the tables
  print(table_sex_vhl)  # Count table
  print(round(prop_sex_vhl * 100, 2))  # Percentage table
  
  # Perform Chi-squared and Fisher's exact tests
  if (all(table_sex_vhl > 0)) {  # Chi-squared test requires non-zero counts
    chisq_test <- chisq.test(table_sex_vhl)
    print(chisq_test)
  } else {
    cat("Chi-squared test not applicable due to zero counts.\n")
  }
  
  fisher_test <- fisher.test(table_sex_vhl)
  print(fisher_test)
  
  # Create Proportion Plot
  proportion_plot <- ggplot(tumors_stage, aes(x = VHLStatus, fill = Sex)) +
    geom_bar(position = "fill") +
    labs(y = "Proportion", title = paste("Stage", stage)) +
    scale_fill_manual(values = c("pink", "lightblue")) +
    theme_minimal()
  proportion_plots[[stage]] <- proportion_plot
  
  # Create Count Plot
  count_plot <- ggplot(tumors_stage, aes(x = VHLStatus, fill = Sex)) +
    geom_bar(position = "stack") +
    labs(y = "Count", title = paste("Stage", stage)) +
    scale_fill_manual(values = c("pink", "lightblue")) +
    theme_minimal()
  count_plots[[stage]] <- count_plot
}

# Display all Proportion Plots in order
if (length(proportion_plots) > 0) {
  cat("\n### Displaying Proportion Plots ###\n")
  gridExtra::grid.arrange(grobs = proportion_plots, ncol = 2)
}

# Display all Count Plots in order
if (length(count_plots) > 0) {
  cat("\n### Displaying Count Plots ###\n")
  gridExtra::grid.arrange(grobs = count_plots, ncol = 2)
}






#### stage ####
# Read the data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')

# Exclude rows with NA values in 'Stage' or 'VHL' columns
tumors_clean <- tumors[tumors$VHL %in% c("1", "0"), ]

# Check the structure of the cleaned dataset
head(tumors_clean)

tumors_clean <- tumors_clean[tumors_clean$VHL == "1", ]

# Now proceed with the rest of the analysis using the cleaned data
# Contingency table for Stage and VHL
table_stage_vhl <- table(tumors_clean$Stage, tumors_clean$KDM5C_SETD2)
table_stage_vhl

# Chi-Square Test
chisq_test_stage <- chisq.test(table_stage_vhl)
print(chisq_test_stage)

# Fisher's Exact Test (if expected frequencies are small)
fisher_test_stage <- fisher.test(table_stage_vhl)
print(fisher_test_stage)

# Convert Stage to an ordered factor (for ordinal data)
tumors_clean$Stage <- factor(tumors_clean$Stage, levels = c("I", "II", "III", "IV", "Missing"), ordered = TRUE)

# Perform Wilcoxon rank-sum test (if stage is ordinal)
wilcox_test_stage <- wilcox.test(as.numeric(Stage) ~ KDM5C_SETD2, data = tumors_clean)
print(wilcox_test_stage)

# Visualization: Proportion bar plot
ggplot(tumors_clean, aes(x = factor(KDM5C_SETD2), fill = Stage)) +  
  geom_bar(position = "fill") +
  labs(x = "KDM5C_SETD2", y = "Proportion", title = "Stage Distribution by VHL Status") +
  scale_x_discrete(labels = c("0" = "BAP1wt", "1" = "BAP1mut")) +  
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()


# Violin plot for ordinal stage comparison
ggplot(tumors_clean, aes(x = PBRM1_PBRM1_2, y = as.numeric(Stage), fill = PBRM1_PBRM1_2)) +
  geom_violin(trim = FALSE) +
  labs(y = "Tumor Stage (Ordinal)", title = "Stage Comparison by VHL Status") +
  scale_fill_manual(values = c("red1", "blue")) +
  theme_minimal()









#### Tumor Size ####
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')
tumors <- tumors[tumors$VHL == "0", ]
#tumors <- tumors[tumors$Stage == "I", ]
# Check the distribution of Tumor Size
tumors_clean <- tumors[!is.na(tumors$KDM5C_SETD2) & !is.na(tumors$TumourSize_raw), ]
shapiro.test(tumors_clean$TumourSize_raw)

# If data is normally distributed, use t-test, else Wilcoxon test
t_test_size <- t.test(TumourSize_raw ~ KDM5C_SETD2, data = tumors_clean)
print(t_test_size)

# Wilcoxon test (if size data is non-normal)
wilcox_test_size <- wilcox.test(TumourSize_raw ~ KDM5C_SETD2, data = tumors_clean)
print(wilcox_test_size)

# Visualization: Boxplot of Tumor Size by PBRM1_SETD2 status
ggplot(tumors_clean, aes(x = factor(KDM5C_SETD2), y = TumourSize_raw, fill = factor(KDM5C_SETD2))) +
  geom_boxplot() +
  labs(y = "Tumor Size (mm)", title = "Tumor Size") +
  scale_fill_manual(values = c("brown2", "royalblue3"), labels = c("Co-occur", "Other")) +
  scale_x_discrete(labels = c("0" = "PBRM1_SETD2wt", "1" = "PBRM1_SETD2mut")) + 
  theme_minimal()




#### Tumor Size Stratified by Stage####
library(ggplot2)
library(gridExtra)
library(readxl)

# Load the data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')


tumors <- tumors[tumors$VHL == "1", ]
# Ensure Stage is a factor and ordered
tumors$Stage <- factor(tumors$Stage, levels = c("I", "II", "III", "IV"))

# Filter out NA values in tumor size and VHL status
tumors_clean <- tumors[!is.na(tumors$BAP1) & !is.na(tumors$TumourSize_raw), ]

# Initialize list to store plots
tumor_size_plots <- list()

# Loop through each stage to create stratified plots
for (stage in levels(tumors_clean$Stage)) {
  cat("\n### Stage:", stage, "###\n")
  
  # Subset data for the current stage
  tumors_stage <- tumors_clean[tumors_clean$Stage == stage, ]
  
  # Check if data exists for the stage
  if (nrow(tumors_stage) == 0) {
    cat("No valid data for Stage:", stage, "\n")
    next
  }
  
  # Perform Normality Test for Tumor Size
  shapiro_test <- shapiro.test(tumors_stage$TumourSize_raw)
  print(shapiro_test)
  
  # Perform t-test or Wilcoxon test based on normality
  if (shapiro_test$p.value > 0.05) {
    t_test_size <- t.test(TumourSize_raw ~ BAP1, data = tumors_stage)
    print(t_test_size)
  } else {
    wilcox_test_size <- wilcox.test(TumourSize_raw ~ BAP1, data = tumors_stage)
    print(wilcox_test_size)
  }
  
  # Create a boxplot for Tumor Size by VHL status
  tumors_stage_filtered <- tumors_stage[tumors_stage$BAP1 %in% c(0, 1), ]
  
  tumor_size_plot <- ggplot(tumors_stage_filtered, aes(x = factor(BAP1), y = TumourSize_raw, fill = factor(BAP1))) +
    geom_boxplot() +
    labs(x = "BAP1 Mutation Status", y = "Tumor Size (mm)", title = paste("Tumor Size - Stage", stage)) +
    scale_x_discrete(labels = c("0" = "BAP1wt", "1" = "BAP1mut")) +  # Modify x-axis labels
    scale_fill_manual(values = c("blue", "red"), labels = c("BAP1wt", "BAP1mut")) +  # Modify legend labels
    theme_minimal()
  

  
  # Store plot in the list
  tumor_size_plots[[stage]] <- tumor_size_plot
}

# Display all Tumor Size Boxplots
if (length(tumor_size_plots) > 0) {
  cat("\n### Displaying Tumor Size Boxplots by Stage ###\n")
  gridExtra::grid.arrange(grobs = tumor_size_plots, ncol = 2)
}




# Create a contingency table
contingency_table <- matrix(c(30, 101, 190, 61), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
fisher_result <- fisher.test(contingency_table)

# Print the results
print(fisher_result)


