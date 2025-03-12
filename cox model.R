#### cox model for Survival ####

# Load required packages
install.packages("survival")
install.packages("survminer")
install.packages("forestmodel") 
library(forestmodel)
library(survival)
library(survminer)

#load file and exclude metastatic sample
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')
tumors <- tumors[tumors$LocalisedRCC == "1", ]
#tumors <- na.omit(tumors)  # delete NA row


# Filtering by VHL status
tumors <- tumors[tumors$VHL == "1", ]
#tumors <- tumors[tumors$Stage == "II", ]

#set reference
tumors$KDM5C_SETD2_4 <- factor(tumors$KDM5C_SETD2_4, levels = c("ZNone", "KDM5C_only", "SETD2_only", "Co-occur"))
tumors$Sex <- factor(tumors$Sex, levels = c("Female", "Male"))
tumors$Age_raw <- as.numeric(tumors$Age_raw)
tumors$Stage <- factor(tumors$Stage, levels = c("I", "II", "III"))

# Define the 5-year DFS
max_time <- 5
tumors$DFSYears_5yr <- pmin(tumors$DFSYears, max_time)  # Cap the DFS time at 5 years
tumors$DFSEvent_5yr <- ifelse(tumors$DFSYears > max_time, 0, tumors$DFSEvent)  # Mark events beyond 5 years as censored

# Create the survival object
surv_object <- Surv(tumors$DFSYears_5yr, tumors$DFSEvent_5yr)

tumors$Co_occurrence <- factor(
  ifelse(tumors$KDM5C == 1 & tumors$SETD2 == 1, "Co-occur",
         ifelse(tumors$KDM5C == 1, "PBRM1_only",
                ifelse(tumors$SETD2 == 1, "SETD2_only", "None"))),
  levels = c("None", "Co-occur", "PBRM1_only", "SETD2_only")  # Ensure correct order
)


# Fit the multivariate Cox Proportional Hazards Model
cox_model <- coxph(surv_object ~ KDM5C_SETD2_4+Age_raw+Sex, data = tumors)

# View model summary
summary(cox_model)

#forest plot
forest_model(cox_model)

#check if the mode is reasonable, if p-value < 0.05, then some variate do not fit
cox.zph(cox_model)







#### Logistic model for Stage #### 
library(nnet)
library(ggplot2)
library(dplyr)
library(broom)

tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')
tumors$Stage <- factor(tumors$Stage, levels = c("I", "II", "III", "IV"))
tumors <- tumors[tumors$VHL == "0", ]
logit_model <- multinom(Stage ~ PBRM1 + SETD2 + KDM5C + BAP1 + TP53, data = tumors)
summary(logit_model)



# save results into data frame
tidy_results <- broom::tidy(logit_model, conf.int = TRUE)

tidy_results <- tidy_results %>%
  filter(term != "(Intercept)")

# Extract coefficients from the multinomial logistic regression model
tidy_results <- broom::tidy(logit_model, conf.int = TRUE)

# Filter out intercepts and ensure Stage II, III, and IV appear as separate rows
tidy_results <- tidy_results %>%
  filter(term != "(Intercept)") %>%
  mutate(y.level = factor(y.level, levels = c("IV", "III", "II")))  # Ensure Stage IV is at the top

# Convert p-values for display
tidy_results <- tidy_results %>%
  mutate(p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))  # Format p-values

# Adjust position for p-value labels
tidy_results <- tidy_results %>%
  mutate(v_adjust = case_when(
    y.level == "II" ~ -0.5,   # Stage II p-value above
    y.level == "III" ~ 1.9,   # Stage III p-value below
    y.level == "IV" ~ 2.0     # Stage IV p-value further below
  ))

# Define colors: Stage II (Orange), Stage III (Blue), Stage IV (Green)
stage_colors <- c("II" = "#E69F00", "III" = "#0072B2", "IV" = "#009E73")

# Plot forest plot with Stage II, III, and IV separately
ggplot(tidy_results, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = y.level)) +
  geom_pointrange(position = position_dodge(width = 0.3), size = 0.6) +  # Separate points for Stage II, III, IV
  geom_text(aes(label = paste0("p=", p_label), vjust = v_adjust), 
            position = position_dodge(width = 0.5), 
            size = 3.5) +  # Adjust p-value placement
  coord_flip() +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Reference line at 0 (log odds)
  labs(title = "Forest Plot of Multinomial Logistic Regression",
       x = "Predictors",
       y = "Log Odds Ratio",
       color = "Stage") +
  scale_color_manual(values = stage_colors)  # Set new colors for all three stages







#### linear model for Tumor Size ####

# Load required libraries
library(ggplot2)
library(dplyr)
library(broom)
library(readxl)

# Load data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')

tumors <- tumors[tumors$VHL == "0", ]
tumors <- tumors[tumors$Stage == "I", ]
# Ensure TumourSize_raw is numeric
tumors$TumourSize_raw <- as.numeric(tumors$TumourSize_raw)

# Set categorical variables as factors
tumors$Sex <- factor(tumors$Sex, levels = c("Female", "Male"))
tumors$Stage <- factor(tumors$Stage, levels = c("I", "II", "III", "IV"))
tumors$PBRM1_SETD2 <- factor(tumors$PBRM1_SETD2, levels = c("Other", "Co-occur"))

tumors$Co_occurrence <- factor(
  ifelse(tumors$PBRM1 == 1 & tumors$SETD2 == 1, "Co-occur",
         ifelse(tumors$PBRM1 == 1, "PBRM1_only",
                ifelse(tumors$SETD2 == 1, "SETD2_only", "None"))),
  levels = c("None", "Co-occur", "PBRM1_only", "SETD2_only")  # Ensure correct order
)

# Optional: Remove missing values
#tumors <- na.omit(tumors)  

# Fit linear regression model predicting Tumor Size
lm_model <- lm(TumourSize_raw ~ PBRM1_SETD2+Sex+Age_raw+BMI, data = tumors)

# View model summary
summary(lm_model)

# Extract coefficients and confidence intervals
tidy_results <- broom::tidy(lm_model, conf.int = TRUE)

# Remove intercept
tidy_results <- tidy_results %>% filter(term != "(Intercept)")

# Convert p-values for display
tidy_results <- tidy_results %>%
  mutate(p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))

# Define colors for significance
tidy_results$significance <- ifelse(tidy_results$p.value < 0.05, "Significant", "Not Significant")

# Ensure term order follows model order
#tidy_results$term <- factor(tidy_results$term, 
#                            levels = c("BAP1", "TP53", "Co_occurrenceSETD2_only", "Co_occurrencePBRM1_only", "Co_occurrenceCo-occur" 
#                                         ))  # Match lm model

# Plot forest plot with correct order
ggplot(tidy_results, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = significance)) +
  geom_pointrange(size = 0.8) +  
  geom_text(aes(label = paste0("p=", p_label)), hjust = -0.1, vjust = -1, size = 3) + # Display p-values
  coord_flip() +  
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Reference line at 0
  labs(title = "Forest Plot of Linear Regression (Tumor Size)",
       x = "Predictors",
       y = "Estimated Effect (β)") +
  scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "#666666"))


forest_model(lm_model)




#### Stage HvsL ####
# Load necessary libraries
library(nnet)
library(ggplot2)
library(dplyr)
library(broom)
library(readxl)

# Load data
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')

tumors <- tumors[tumors$VHL == "0", ]

# Check for NA values and remove them if needed
#tumors <- tumors %>% filter(!is.na(StageHvsL))

# Convert to factor and explicitly set levels
tumors$StageHvsL <- factor(tumors$StageHvsL, levels = c("Low", "High"))

# Verify counts again
table(tumors$StageHvsL)

# Fit multinomial logistic regression
logit_model <- multinom(StageHvsL ~ PBRM1 + SETD2 + KDM5C + BAP1 + TP53, data = tumors)

summary(logit_model)

# Extract coefficients and confidence intervals
tidy_results <- broom::tidy(logit_model, conf.int = TRUE)

# Remove intercept term
tidy_results <- tidy_results %>%
  filter(term != "(Intercept)")

# Convert p-values for display
tidy_results <- tidy_results %>%
  mutate(p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))

# Create a significance column
tidy_results <- tidy_results %>%
  mutate(significance = ifelse(p.value < 0.05, "Significant", "Not Significant"))

# Adjust vertical position of p-values
tidy_results <- tidy_results %>%
  mutate(v_adjust = -0.5)  # Shift p-value labels downward

# Plot forest plot
ggplot(tidy_results, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = significance)) +
  geom_pointrange(position = position_dodge(width = 0.3), size = 0.6) +  # Plot log-odds ratio
  geom_text(aes(label = paste0("p=", p_label)), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5, size = 3.5) +  # Adjust p-value position
  coord_flip() +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Reference line at 0 (log odds)
  labs(title = "Forest Plot of Multinomial Logistic Regression (StageHvsL)",
       x = "Predictors",
       y = "Log Odds Ratio") +  
  scale_color_manual(values = c("Significant" = "purple", "Not Significant" = "gray")) 

