library(survival)
library(survminer)
library(readxl)

# Read in the data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')
#exclude metastatic samples
tumors <- tumors[tumors$LocalisedRCC == "1", ]

# Filtering
tumors <- tumors[tumors$VHL == "1", ]
tumors <- tumors[tumors$Stage == "II", ]

# Define the 5-year DFS
max_time <- 5
tumors$DFSYears_5yr <- pmin(tumors$DFSYears, max_time)  # Cap the DFS time at 5 years
tumors$DFSEvent_5yr <- ifelse(tumors$DFSYears > max_time, 0, tumors$DFSEvent)  # Mark events beyond 5 years as censored

# Define the survival object
tumorSurv_5yr <- Surv(time = tumors$DFSYears_5yr, event = tumors$DFSEvent_5yr)

# Fit the Kaplan-Meier survival model
tumorKM_5yr <- survfit(tumorSurv_5yr ~ KDM5C_SETD2_4, data = tumors, type = "kaplan-meier")

# Perform a log-rank test
logrank_test_5yr <- survdiff(tumorSurv_5yr ~ KDM5C_SETD2_4, data = tumors)

# Calculate the p-value for the log-rank test
logrank_pval_5yr <- 1 - pchisq(logrank_test_5yr$chisq, df = length(logrank_test_5yr$n) - 1)

# Print the log-rank p-value for reference
print(paste("5-year log-rank p-value:", round(logrank_pval_5yr, 3)))


ggsurvplot(
  tumorKM_5yr,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval_5yr, 3)),
  pval.coord = c(3, 1),  # Add the p-value to the plot
  risk.table = TRUE,     # Enable the risk table at the bottom
  legend.labs = c("0", "1", "2", "3","4","5","6"),
  legend = c(0.15, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("purple", "olivedrab3", "steelblue1", "azure4","gold","red","blue"),
  xlab = "Survival (years)",
  ylab = "Proportion surviving"
)


ggsurvplot(
  tumorKM_5yr,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval_5yr, 3)),
  pval.coord = c(2.8, 0.8),  # Add the p-value to the plot
  risk.table = TRUE,     # Enable the risk table at the bottom
  legend.labs = c("KDM5C only", "SETD2 only", "None"),
  legend = c(0.15, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("olivedrab3", "steelblue1", "azure4"),
  xlab = "Survival (years)",
  ylab = "Proportion surviving"
)


# Plot the Kaplan-Meier curve with the p-value and risk table
# 4 curves: Both, Gene1, Gene2, None
ggsurvplot(
  tumorKM_5yr,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval_5yr, 3)),
  pval.coord = c(3, 1),  # Add the p-value to the plot
  risk.table = TRUE,     # Enable the risk table at the bottom
  legend.labs = c("KDM5C_SETD2", "KDM5C only", "SETD2 only", "None"),
  legend = c(0.15, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("purple", "olivedrab3", "steelblue1", "azure4"),
  xlab = "Survival (years)",
  ylab = "Proportion surviving"
)


# 2 curves: Both Genes, Other Circunstances
ggsurvplot(
  tumorKM_5yr,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval_5yr, 3)),
  pval.coord = c(4, 0.98),# Add the p-value to the plot
  risk.table = TRUE, # Enable the risk table at the bottom
  legend.labs = c("KDM5C_SETD2_4", "Other"),
  legend = c(0.13, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("brown2", "royalblue3"),
  xlab = "Survival (years)",
  ylab = "Disease-Free Survival"
)

