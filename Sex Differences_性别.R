#### Sex Differences ####

# load input data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')

# Summarize counts by Sex and VHL Status
table_sex_vhl <- table(tumors$Sex, tumors$VHL)

# Calculate proportions
prop_sex_vhl <- prop.table(table_sex_vhl, margin = 2)

# Print the tables
print(table_sex_vhl)  # Count table
print(round(prop_sex_vhl * 100, 2))  # Percentage table

# Pearson's Chi-squared test with Yates' continuity correction
chisq_test <- chisq.test(table_sex_vhl)
print(chisq_test)

# Fisher's Exact Test for Count Data
fisher_test <- fisher.test(table_sex_vhl)
print(fisher_test)

library(ggplot2)

# Proportion
ggplot(tumors, aes(x = VHL, fill = Sex)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", title = "Sex Proportion by VHL Status") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal()

# Counts
ggplot(tumors, aes(x = VHL, fill = Sex)) +
  geom_bar(position = "stack") +
  labs(y = "Count", title = "Sex Counts by VHL Status") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal()

