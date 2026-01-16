
#All Packages
library(readxl)
library(dplyr) 
library(tidyverse) 
library(writexl)
library(dplyr)
library(missForest)


#Read all input files
country_income_status <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/Country_Income_Status_WHO.xlsx', sheet = "List of economies")

rep_wash_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/rep_wash/wash_data.xlsx', sheet = "Cleaned Data")
env_health_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/rep_gho_env/env_health_data.xlsx', sheet = "cleaned_data")
health_spending_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/HealthSpending/healthspeinding_cleaned.xlsx', sheet = "Sheet 1")
healthcare_access_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/healthcare_access_rep_dhs_hca/healthcare_access_data.xlsx', sheet = "Cleaned Data")
health_determinant_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/Health_Determinants_rep_dhs_unicef_sdh/healthdeterminants_cleaned.xlsx', sheet = "Sheet 1")
climate_cleaned_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/climate_cleaned_data.xlsx', sheet = "Sheet 1")
development_indices_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/development_indices_globalDataLab_rep_gdl1/development_indices_subset.xlsx', sheet = "CleanedData")
antibioticconsumption_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/antibioticconsumption_consumption.xlsx', sheet = "Sheet 1")
adult_health_nutrition <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/adult_health_nutrition_rep_dhs_ahn/adult_health_nutrition_data.xlsx', sheet = "cleaned data")

##Global Governance
gg_ControlofCorruption_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/global_governance_cleaned_data/ControlofCorruption.xlsx', sheet = "Sheet 1")
gg_GovernmentEffectiveness_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/global_governance_cleaned_data/GovernmentEffectiveness.xlsx', sheet = "Sheet 1")
gg_PoliticalStabilityNoViolence_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/global_governance_cleaned_data/PoliticalStabilityNoViolence.xlsx', sheet = "Sheet 1")
gg_RegulatoryQuality_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/global_governance_cleaned_data/RegulatoryQuality.xlsx', sheet = "Sheet 1")
gg_RuleofLaw_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/global_governance_cleaned_data/RuleofLaw.xlsx', sheet = "Sheet 1")
gg_VoiceandAccountability_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/global_governance_cleaned_data/VoiceandAccountability.xlsx', sheet = "Sheet 1")


#WHO Data
who_data <- read_excel('/Users/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/whocleaned_data.xlsx', sheet = "Sheet 1")

# Create a list of your data frames
data_frames_list <- list(gg_VoiceandAccountability_data, 
                         gg_RuleofLaw_data, 
                         gg_RegulatoryQuality_data, 
                         gg_PoliticalStabilityNoViolence_data, 
                         gg_GovernmentEffectiveness_data,
                         gg_ControlofCorruption_data, 
                         adult_health_nutrition, 
                         antibioticconsumption_data, 
                         development_indices_data, 
                         climate_cleaned_data, 
                         health_determinant_data, 
                         healthcare_access_data, 
                         health_spending_data, 
                         env_health_data, 
                         rep_wash_data)



# Ensure all data frames have the correct data types
data_frames_list_corrected <- lapply(data_frames_list, function(df) {
  df %>%
    mutate(across(c(source, iso3, indicator), as.character), # Convert to character
           across(c(year, value), as.numeric)) # Convert to numeric (double)
})

# Combine all corrected data frames into one
combined_data_frame <- bind_rows(data_frames_list_corrected)

# Inspect the first few rows of the combined data frame
head(combined_data_frame)
nrow(combined_data_frame)

# Extract distinct values
variables_with_sources <- combined_data_frame %>%
  distinct(indicator, source)

# Save to an Excel file
write_xlsx(variables_with_sources, "/Users/nikhilg/Documents/Vivli\ Research\ Project/R\ Output/distinct_indicator_source.xlsx")


# Drop the "source" column from the dataframe
combined_data_frame <- combined_data_frame %>%
  select(-source)

print(combined_data_frame)

# Pivot the combined data frame with aggregation on average
wide_data_frame <- pivot_wider(combined_data_frame, 
                               names_from = indicator, 
                               values_from = value,
                               values_fn = list(value = ~mean(.x, na.rm = TRUE))) # Aggregate using average, ignore NA

# Inspect the first few rows of the widened data frame
head(wide_data_frame)

# Adjust data types in who_data
who_data <- who_data %>%
  mutate(iso3 = as.character(iso3),
         year = as.double(year))

# Adjust data types in wide_data_frame
wide_data_frame <- wide_data_frame %>%
  mutate(iso3 = as.character(iso3),
         year = as.double(year))

# Combine the dataframes based on 'iso3' and 'year'
# Perform full join and keep only columns from who_data when duplicates exist
all_input_variables <- full_join(
  who_data, 
  wide_data_frame, 
  by = c("iso3", "year"),
  suffix = c("", "_drop") # Add suffix to conflicting columns from wide_data_frame
) %>%
  select(-ends_with("_drop"))  # Remove duplicate columns from wide_data_frame

# View the first few rows of the combined dataframe
head(all_input_variables)

# Select only the key and the columns you want from country_income_status
country_income_selected <- country_income_status[, c("Code", "Income group", "Region")]

# Perform the left join
result_dataframe <- merge(all_input_variables, country_income_selected, 
                          by.x = "iso3", by.y = "Code", 
                          all.x = TRUE)


# Replace all NaN values with NA in the entire DataFrame
# Iterate over each column of the DataFrame
for (col_name in names(result_dataframe)) {
  # Apply the replacement only to numeric columns
  if (is.numeric(result_dataframe[[col_name]])) {
    result_dataframe[[col_name]][is.nan(result_dataframe[[col_name]])] <- NA
  }
}

# Remove rows where income_group is blank (either "" or NA)
result_dataframe_cleaned <- result_dataframe %>%
  filter(`Income group` != "" & !is.na(`Income group`))

# View the cleaned dataframe
nrow(result_dataframe_cleaned)

# Sorting by iso3 and then by year, both in ascending order
result_dataframe_sorted <- result_dataframe_cleaned %>%
  arrange(iso3, year)

#Only select data where year > 2004
result_dataframe_sorted <- subset(result_dataframe_sorted, year > 2004)

print(result_dataframe_sorted)


# Apply missForest imputation separately for 'iso3', then 'Region', then 'Income group'
result_dataframe_sorted_test <- result_dataframe_sorted %>%
  group_by(iso3, Region, `Income group`) %>%  # Group by multiple levels
  group_modify(~ {
    df <- as.data.frame(.x)  # Convert to data frame (required for missForest)
    
    # Identify numeric columns
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    
    # Subset only numeric columns for imputation
    df_numeric <- df[numeric_cols]
    
    # Apply missForest only if there's at least one missing value in numeric columns
    if (any(is.na(df_numeric))) {
      df_numeric_imputed <- missForest(df_numeric)$ximp  # Impute missing values
    } else {
      df_numeric_imputed <- df_numeric  # No missing values, return original numeric columns
    }
    
    # Replace the numeric columns in the original dataset with the imputed values
    df[numeric_cols] <- df_numeric_imputed
    
    return(df)
  }) %>%
  ungroup()



# Calculate the percentage of missing values for each column
missing_percentages <- colSums(is.na(result_dataframe_sorted_test)) / nrow(result_dataframe_sorted_test) * 100

# Sort the percentages in descending order
missing_percentages_sorted <- sort(missing_percentages, decreasing = TRUE)

# Print the sorted percentages
print(missing_percentages_sorted)

#Remove the column named "Not Required"
result_dataframe_sorted_test <- result_dataframe_sorted_test %>%
  select(-`Not Required`)

# Save final input indicator data file
write_xlsx(result_dataframe_sorted_test, "/Users/nikhilg/Documents/Vivli\ Research\ Project/R\ Output/cleaned_variable_dataset.xlsx")


