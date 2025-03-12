#### Tumor Stages ####


# Exclude rows with NA values in 'Stage' or 'VHL' columns
tumors_clean <- tumors[!is.na(tumors$Stage) & !is.na(tumors$VHL), ]

# Check the structure of the cleaned dataset
head(tumors_clean)

# Contingency table for Stage and VHL
table_stage_vhl <- table(tumors_clean$Stage, tumors_clean$VHL)
table_stage_vhl

# Convert Stage to an ordered factor (for ordinal data)
tumors_clean$Stage <- factor(tumors_clean$Stage, levels = c("I", "II", "III", "IV"), ordered = TRUE)

# Chi-Square Test
chisq_test_stage <- chisq.test(table_stage_vhl)
print(chisq_test_stage)

# Fisher's Exact Test (if expected frequencies are small)
fisher_test_stage <- fisher.test(table_stage_vhl)
print(fisher_test_stage)

# Perform Wilcoxon rank-sum test (if stage is ordinal)
wilcox_test_stage <- wilcox.test(as.numeric(Stage) ~ VHL, data = tumors_clean)
print(wilcox_test_stage)


# Proportion bar plot
library(ggplot2)
ggplot(tumors_clean, aes(x = VHL, fill = Stage)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", title = "Stage Distribution by VHL Status") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()

# Count bar plot
library(ggplot2)
ggplot(tumors_clean, aes(x = VHL, fill = Stage)) +
  geom_bar(position = "stack") +
  labs(y = "Proportion", title = "Stage Distribution by VHL Status") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()

