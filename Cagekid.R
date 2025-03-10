#### Sex Differnce ####
tumors <- read_excel('/Users/pliu/Desktop/Cagekid/Cagekid_MAF_ALL.xlsx', sheet = 'Cagekid_clin')
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


#==========================================================================
#### Stages ####
# Exclude rows with NA values in 'Stage' or 'VHL' columns
tumors_clean <- tumors[!is.na(tumors$Stage) & !is.na(tumors$VHL), ]

# Check the structure of the cleaned dataset
head(tumors_clean)

# Now proceed with the rest of the analysis using the cleaned data
# Contingency table for Stage and VHL
table_stage_vhl <- table(tumors_clean$Stage, tumors_clean$VHL)
table_stage_vhl

# Chi-Square Test
chisq_test_stage <- chisq.test(table_stage_vhl)
print(chisq_test_stage)

# Fisher's Exact Test (if expected frequencies are small)
fisher_test_stage <- fisher.test(table_stage_vhl)
print(fisher_test_stage)

# Convert Stage to an ordered factor (for ordinal data)
tumors_clean$Stage <- factor(tumors_clean$Stage, levels = c("I", "II", "III", "IV"), ordered = TRUE)

# Perform Wilcoxon rank-sum test (if stage is ordinal)
wilcox_test_stage <- wilcox.test(as.numeric(Stage) ~ VHL, data = tumors_clean)
print(wilcox_test_stage)

# Visualization: Proportion bar plot
library(ggplot2)
ggplot(tumors_clean, aes(x = VHL, fill = Stage)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", title = "Stage Distribution by VHL Status") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()

# Violin plot for ordinal stage comparison
ggplot(tumors_clean, aes(x = VHL, y = as.numeric(Stage), fill = VHL)) +
  geom_violin(trim = FALSE) +
  labs(y = "Tumor Stage (Ordinal)", title = "Stage Comparison by VHL Status") +
  scale_fill_manual(values = c("red1", "blue")) +
  theme_minimal()

# Check if any other analysis needs to be done after cleaning the data



#============================================================================
#### Tumor Size ####
tumors <- read_excel('/Users/jessie/Desktop/Cagekid_MAF_ALL.xlsx', sheet = 'Cagekid_clin')
tumors <- tumors[tumors$Stage == "IV", ]
tumors <- tumors[tumors$VHL == "VHLmut", ]
tumors <- tumors[tumors$Sex == "Male", ]
# Check the distribution of Tumor Size
tumors_clean <- tumors[!is.na(tumors$VHL) & !is.na(tumors$TumourSize_raw), ]
shapiro.test(tumors_clean$TumourSize_raw)

# If data is normally distributed, use t-test, else Wilcoxon test
t_test_size <- t.test(TumourSize_raw ~ Sex, data = tumors_clean)
print(t_test_size)

# Wilcoxon test (if size data is non-normal)
wilcox_test_size <- wilcox.test(TumourSize_raw ~ Sex, data = tumors_clean)
print(wilcox_test_size)

# Visualization: Boxplot of Tumor Size by VHL status
ggplot(tumors_clean, aes(x = Sex, y = TumourSize_raw, fill = Sex)) +
  geom_boxplot() +
  labs(y = "Tumor Size (mm)", title = "Tumor Size") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal()







# Initialize a list to store results
results <- list()

# Initialize a vector to store p-values
p_values <- c()

# Loop through each Stage group
for (group in levels(tumors_clean$Stage)) {
  
  # Subset the data for the current Stage group
  subset_data <- tumors_clean[tumors_clean$Stage == group, ]
  
  # Check the distribution of Tumor Size (Shapiro-Wilk test)
  shapiro_test <- shapiro.test(subset_data$TumourSize_raw)
  
  # Perform the appropriate test based on normality
  if (shapiro_test$p.value > 0.05) {
    # Normal distribution: t-test
    test_result <- t.test(TumourSize_raw ~ Sex, data = subset_data)
  } else {
    # Non-normal distribution: Wilcoxon test
    test_result <- wilcox.test(TumourSize_raw ~ Sex, data = subset_data)
  }
  
  # Store the result for the current Stage group
  results[[group]] <- list(
    "Shapiro-Wilk p-value" = shapiro_test$p.value,
    "Test p-value" = test_result$p.value,
    "Test Method" = ifelse(shapiro_test$p.value > 0.05, "t-test", "Wilcoxon test")
  )
  
  # Add the p-value to the vector
  p_values <- c(p_values, test_result$p.value)
}

# Print all p-values
print(p_values)


library(dplyr)
library(ggplot2)

# Calculate the counts for each combination of Stage group and KDM5C_PBRM1
count_data <- tumors_clean %>%
  group_by(Stage, Sex) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(group_number = as.integer(factor(Stage)))  # Add Stage group number

# Calculate Wilcoxon test p-values for each Stage group
p_values_data <- tumors_clean %>%
  group_by(Stage) %>%
  summarise(
    p_value = wilcox.test(TumourSize_raw ~ Sex, data = .)$p.value,
    .groups = "drop"
  )

# Combine count_data and p_values_data for plotting
annotation_data <- count_data %>%
  left_join(p_values_data, by = "Stage")

# Plot the distribution of Tumor Size by v and Stage group
ggplot(tumors_clean, aes(x = Sex, y = TumourSize_raw, fill = Sex)) +
  geom_boxplot() +
  facet_wrap(~ Stage, scales = "free_y") +  # Facet by Stage groups
  labs(y = "Tumor Size (mm)", title = "Tumor Size Comparison by Stage Group") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal() +
  # Add sample counts to the plot
  geom_text(data = annotation_data, aes(x = Sex, 
                                        y = max(tumors_clean$TumourSize_raw) * 1.1, 
                                        label = paste("n =", count)),
            position = position_dodge(width = 0.75), size = 3, vjust = 0) +
  # Add Wilcoxon p-values to the plot
  geom_text(data = p_values_data, aes(x = 1.5, 
                                      y = max(tumors_clean$TumourSize_raw) * 1.2, 
                                      label = paste("p =", signif(p_value, 3))),
            inherit.aes = FALSE, size = 3, vjust = 0)

print(p_values_data)


#============================================================================
#### Survival Analysis ####
install.packages("survminer")
# Load required libraries
library(survival)
library(survminer)
library(readxl)

# Read in the data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid_MAF_ALL.xlsx', sheet = 'Cagekid_clin')

# Check the structure of the data
head(tumors)
tumors <- tumors[tumors$VHL == "VHLwt", ]
#tumors <- tumors[tumors$Sex == "Male", ]
# Define the survival object
tumorSurv <- Surv(time = tumors$DFS_5Years, event = tumors$Dead)

# Fit the Kaplan-Meier survival model
tumorKM <- survfit(tumorSurv ~ KDM5C, data = tumors, type = "kaplan-meier")

# Perform a log-rank test
logrank_test <- survdiff(tumorSurv ~ KDM5C, data = tumors)

# Calculate the p-value for the log-rank test
logrank_pval <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)

# Print the log-rank p-value for reference
print(paste("Log-rank p-value:", round(logrank_pval, 3)))

# Plot the Kaplan-Meier curve with the p-value and risk table
ggsurvplot(
  tumorKM,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval, 3)),
  pval.coord = c(3, 0.95),# Add the p-value to the plot
  risk.table = TRUE,                            # Enable the risk table at the bottom
  legend.labs = c("Other", "KDM5C"),
  legend = c(0.8, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("gold", "purple"),
  xlab = "Survival (years)",
  ylab = "Proportion surviving"
)


ggsurvplot(
  tumorKM,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval, 3)),
  pval.coord = c(3, 1),# Add the p-value to the plot
  risk.table = TRUE,                            # Enable the risk table at the bottom
  legend.labs = c("KDM5C_PBRM1", "KDM5C only","PBRM1 only", "None"),
  legend = c(0.2, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("purple","olivedrab3","steelblue1","azure4"),
  xlab = "Survival (years)",
  ylab = "Proportion surviving"
)


#============================================================================
#### Fisher's exact test ####
# Create the contingency table
data <- matrix(c(454, 127, 268, 91), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(data)

# Print the results
print(result)



#============================================================================
#### Stage ####
tumors <- read_excel('/Users/pliu/Desktop/Cagekid/Cagekid_MAF_ALL.xlsx', sheet = 'Cagekid_clin')
tumors <- tumors[tumors$VHL == "VHLmut", ]
tumors <- tumors[!is.na(tumors$Stage) & !is.na(tumors$TumourSize_raw), ]
shapiro.test(tumors_clean$TumourSize_raw)
library(ggplot2)
ggplot(tumors, aes(x = Sex, fill = as.factor(Stage))) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", title = "Mutation Status vs Tumor Stage") +
  scale_x_discrete(labels = c("Both" = "Both", "PBRM1" = "PBRM1", "SETD2" = "SETD2", "YYY" = "Neither")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

tumors$Stage <- factor(tumors$Stage, levels = c("I", "II", "III", "IV"), ordered = TRUE)
ggplot(tumors, aes(x = Sex, y = as.numeric(Stage), fill = Sex)) +
  geom_violin(trim = FALSE) +
  labs(y = "Tumor Stage (Ordinal)", title = "Stage Comparison by Sex") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  scale_y_continuous(breaks = seq(0, max(as.numeric(tumors$Stage)), by = 1)) +  # Set y-axis breaks
  theme_minimal()




#### Ordinal Logistic Regression ####
# Load necessary libraries
library(MASS)        # For ordinal logistic regression
library(readxl)  # For fast MAF file reading
library(dplyr)       # For data manipulation

# Step 1: Load the MAF file
# Replace "your_file.maf" with the actual path to your MAF file
tumors <- read_excel('/Users/jessie/Desktop/Cagekid_MAF_ALL.xlsx', sheet = 'Cagekid_clin')

# Step 2: Preprocess the MAF data
# Ensure the file has the necessary columns for Stage and VHL mutation status
# Modify 'Tumor_Stage' and 'VHL_Status' to match the actual column names in your MAF file
# Example: Tumor_Stage = "I, II, III, IV"; VHL_Status = "mutated" or "wildtype"
tumors <- tumors %>% filter(!is.na(Stage), !is.na(VHL))


# Step 3: Prepare data
# Convert Tumor_Stage to an ordered factor
tumors$Stage <- factor(tumors$Stage, levels = c("I", "II", "III", "IV"), ordered = TRUE)

# Ensure VHL_Status is a factor (categorical variable)
tumors$VHL <- factor(tumors$VHL, levels = c("VHLwt", "VHLmut"))

# Step 4: Fit ordinal logistic regression model
model <- polr(Stage ~ VHL, data = tumors, Hess = TRUE)

# Step 5: Summary and p-values
summary(model)

# Extract coefficients and calculate p-values
coefficients <- coef(summary(model))
p_values <- pnorm(abs(coefficients[, "t value"]), lower.tail = FALSE) * 2

# Add p-values to coefficients table
coefficients <- cbind(coefficients, "p value" = p_values)
print(coefficients)

# Step 6: Odds ratios (OR) and confidence intervals
OR <- exp(coef(model))
CI <- exp(confint(model))
OR_CI <- cbind(OR, CI)
print(OR_CI)

# Optional: Visualize Tumor Stage Distribution by VHL Status
library(ggplot2)
ggplot(tumors, aes(x = VHL, fill = Stage)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", title = "Tumor Stage Distribution by VHL Status") +
  theme_minimal()
