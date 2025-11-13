# =============================================================================
# INTRODUCTORY PROBABILITY AND STATISTICS: BIOMARKERS AND PAIN ANALYSIS
# =============================================================================

# Load required packages
library(tidyverse)
library(readxl)
library(broom)
library(knitr)

# Set seed for reproducibility
set.seed(123)

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Read the data files
biomarkers_data <- read_excel("C:/Users/swyam/OneDrive/Desktop/Prob and Stats/biomarkers.xlsx")
covariates_data <- read_excel("C:/Users/swyam/OneDrive/Desktop/Prob and Stats/covariates.xlsx")

# Clean and merge the data
merged_data <- biomarkers_data %>%
  separate(Biomarker, into = c("PatientID", "Timepoint"), sep = "-") %>%
  mutate(PatientID = as.numeric(PatientID)) %>%
  left_join(covariates_data, by = "PatientID")

# Clean column names for easier handling
cleaned_data <- merged_data %>%
  rename(
    Sex = "Sex (1=male, 2=female)",
    Smoker = "Smoker (1=yes, 2=no)",
    VAS_inclusion = "VAS-at-inclusion",
    VAS_12months = "Vas-12months"
  )

# =============================================================================
# TASK 1: STATISTICAL HYPOTHESIS TESTING
# =============================================================================

# Research Question: Do biomarker levels at inclusion vary between males and females?

# Filter for baseline data only (0weeks)
baseline_data <- cleaned_data %>%
  filter(Timepoint == "0weeks")

# Define biomarker names
biomarkers <- c("IL-8", "VEGF-A", "OPG", "TGF-beta-1", "IL-6", "CXCL9", "CXCL1", "IL-18", "CSF-1")

# Function to perform hypothesis test for a single biomarker
test_biomarker <- function(biomarker_name, data) {
  
  # Check normality using Shapiro-Wilk test
  male_data <- data %>% filter(Sex == 1) %>% pull(biomarker_name)
  female_data <- data %>% filter(Sex == 2) %>% pull(biomarker_name)
  
  shapiro_male <- shapiro.test(male_data)
  shapiro_female <- shapiro.test(female_data)
  
  # Use t-test if both groups are normal, otherwise use Mann-Whitney test
  if (shapiro_male$p.value > 0.05 & shapiro_female$p.value > 0.05) {
    test_result <- t.test(get(biomarker_name) ~ Sex, data = data)
    test_type <- "t-test"
  } else {
    test_result <- wilcox.test(get(biomarker_name) ~ Sex, data = data)
    test_type <- "Mann-Whitney"
  }
  
  # Extract results
  result <- tidy(test_result) %>%
    mutate(
      Biomarker = biomarker_name,
      Test_Type = test_type,
      Shapiro_Male_p = shapiro_male$p.value,
      Shapiro_Female_p = shapiro_female$p.value
    )
  
  return(result)
}

# Perform tests for all biomarkers
results_list <- map(biomarkers, ~test_biomarker(.x, baseline_data))
hypothesis_results <- bind_rows(results_list)

# Apply Bonferroni correction
alpha <- 0.05
n_tests <- length(biomarkers)
bonferroni_alpha <- alpha / n_tests

hypothesis_results <- hypothesis_results %>%
  mutate(
    Significant_Uncorrected = p.value < alpha,
    Significant_Bonferroni = p.value < bonferroni_alpha,
    Bonferroni_Threshold = bonferroni_alpha
  )

# Display results in a clean table
hypothesis_results %>%
  select(Biomarker, Test_Type, statistic, p.value, Significant_Uncorrected, Significant_Bonferroni) %>%
  kable(digits = 4)

# =============================================================================
# TASK 2: REGRESSION MODELLING
# =============================================================================

# Prepare data for regression (baseline measurements only)
regression_data <- cleaned_data %>%
  filter(Timepoint == "0weeks") %>%
  drop_na(VAS_12months)  # Remove rows with missing outcome

# Split data into training (80%) and test (20%) sets
sample_size <- floor(0.8 * nrow(regression_data))
train_indices <- sample(seq_len(nrow(regression_data)), size = sample_size)

train_data <- regression_data[train_indices, ]
test_data <- regression_data[-train_indices, ]

# Build regression model using all biomarkers and covariates
model <- lm(VAS_12months ~ 
              `IL-8` + `VEGF-A` + OPG + `TGF-beta-1` + `IL-6` + 
              CXCL9 + CXCL1 + `IL-18` + `CSF-1` +
              Age + factor(Sex) + factor(Smoker) + VAS_inclusion,
            data = train_data)

# Display model summary
model_summary <- summary(model)
model_summary

# Create a nice table of coefficients
coefficients_table <- tidy(model) %>%
  mutate(
    significant = p.value < 0.05,
    p.value = round(p.value, 4)
  )

coefficients_table %>% kable(digits = 4)

# Model fit statistics
cat("Training Set R-squared:", round(model_summary$r.squared, 4), "\n")
cat("Adjusted R-squared:", round(model_summary$adj.r.squared, 4), "\n")

# =============================================================================
# MODEL EVALUATION ON TEST SET
# =============================================================================

# Make predictions on test set
predictions <- predict(model, newdata = test_data)

# Create evaluation dataframe
evaluation_df <- data.frame(
  Actual = test_data$VAS_12months,
  Predicted = predictions
)

# Calculate performance metrics
mae <- mean(abs(evaluation_df$Actual - evaluation_df$Predicted))
rmse <- sqrt(mean((evaluation_df$Actual - evaluation_df$Predicted)^2))
r_squared_test <- cor(evaluation_df$Actual, evaluation_df$Predicted)^2

cat("\nOUT-OF-SAMPLE EVALUATION METRICS:\n")
cat("Mean Absolute Error (MAE):", round(mae, 4), "\n")
cat("Root Mean Squared Error (RMSE):", round(rmse, 4), "\n")
cat("R-squared on Test Set:", round(r_squared_test, 4), "\n")

# Create prediction plot
prediction_plot <- ggplot(evaluation_df, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual VAS-12months",
       subtitle = "Out-of-sample evaluation",
       x = "Actual VAS-12months",
       y = "Predicted VAS-12months") +
  theme_minimal()

print(prediction_plot)

# =============================================================================
# MULTIPLE TESTING PROBABILITY CALCULATION (for Task 1d)
# =============================================================================

# Probability of at least one Type I error with 9 independent tests
prob_type1_error <- 1 - (1 - 0.05)^9

cat("\nMULTIPLE TESTING ANALYSIS:\n")
cat("Number of tests:", n_tests, "\n")
cat("Probability of at least one Type I error (if all null hypotheses are true):", 
    round(prob_type1_error, 4), "\n")
cat("Bonferroni corrected alpha:", round(bonferroni_alpha, 4), "\n")

# =============================================================================
# SAVE RESULTS FOR REPORT
# =============================================================================

# Save key results to CSV files for easy inclusion in report
write_csv(hypothesis_results, "hypothesis_testing_results.csv")
write_csv(coefficients_table, "regression_coefficients.csv")
write_csv(evaluation_df, "prediction_evaluation.csv")

cat("\nAnalysis complete! Results saved to CSV files.\n")
