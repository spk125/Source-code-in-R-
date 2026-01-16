#All Packages
library(readxl)
library(dplyr) 
library(tidyverse) 
library(writexl)
library(dplyr) 
library(caret)
library(gbm)
library(randomForest)
library(Cubist)
library(xgboost)
library(mgcv)
library(pscl)
library(Metrics)
library(DiceKriging)   # Gaussian Process Regression
library(cluster)  # For spatial clustering
library(groupdata2)  # For group k-fold splitting
library(lme4)
library(car)
library(boot)  # Load bootstrapping library



#Read all AMR Input Data
abaumanii_aminoglycosides_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Abaumanii_dp_Aminoglycosides_final.xlsx', sheet = "Sheet 1")
abaumanii_carbapenems_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Abaumanii_dp_Carbapenems_final.xlsx', sheet = "Sheet 1")
abaumanii_fluoroquinolones_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Abaumanii_dp_Fluoroquinolones_final.xlsx', sheet = "Sheet 1")
escherichiacoli_carbapenems_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Escherichiacoli_dp_Carbapenems_final.xlsx', sheet = "Sheet 1")
escherichiacoli_cephalosporins_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Escherichiacoli_dp_Third_generation_cephalosporins_final.xlsx', sheet = "Sheet 1")
klebsiellapneumoniae_carbapenems_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Klebsiellapneumoniae_dp_Carbapenems_final.xlsx', sheet = "Sheet 1")
klebsiellapneumoniae_cephalosporins_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Klebsiellapneumoniae_dp_Third_generation_cephalosporins_final.xlsx', sheet = "Sheet 1")
pseudomonasaeroginosa_aminoglycosides_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Pseudomonasaeroginosa_dp_Aminoglycosides_final.xlsx', sheet = "Sheet 1")
pseudomonasaeroginosa_cephalosporins_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ amr\ data/Pseudomonasaeroginosa_dp_Third_generation_cephalosporins_final.xlsx', sheet = "Sheet 1")


print(abaumanii_aminoglycosides_data)

#Read indicator data
input_indicator_data <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/R\ Output/cleaned_variable_dataset.xlsx', sheet = "Sheet1")

#Read country longitude, latitude data
country_dim <- read_excel('/Users/amr_user/Documents/Vivli\ Research\ Project/vivli\ indicator\ input\ data/CountryDimData.xlsx', sheet = "Sheet1")
country_dim

# List of data-frame names
df_names <- c("abaumanii_aminoglycosides_data", "abaumanii_carbapenems_data", 
              "abaumanii_fluoroquinolones_data", "escherichiacoli_carbapenems_data", 
              "escherichiacoli_cephalosporins_data", "klebsiellapneumoniae_carbapenems_data", 
              "klebsiellapneumoniae_cephalosporins_data", "pseudomonasaeroginosa_aminoglycosides_data", 
              "pseudomonasaeroginosa_cephalosporins_data")

# New column names
new_col_names <- c("country", "iso3", "year", "TotalResistant", 
                   "TotalSusceptible", "Total", "Resistant_per")

# Loop through each dataframe and rename columns
for (df in df_names) {
  if (exists(df)) {
    assign(df, setNames(get(df), new_col_names))
  } else {
    warning(paste("Dataframe", df, "does not exist in the environment."))
  }
}


# Loop through each dataframe and merge with input_indicator_data
for (df in df_names) {
  if (exists(df)) {
    # Perform merge
    merged_df <- merge(input_indicator_data, get(df)[, c("iso3", "year", "Resistant_per", "Total")], 
                       by = c("iso3", "year"), all.x = TRUE)
    merged_df  <- merge(merged_df , country_dim[, c("iso3", "Latitude", "Longitude")], by = "iso3", all.x = TRUE)
    merged_df <- merged_df[!is.na(merged_df$Latitude) & !is.na(merged_df$Longitude), ]
    merged_df <- merged_df[, c("iso3", "year", "Latitude", "Longitude", setdiff(names(merged_df), c("iso3", "year", "Latitude", "Longitude")))]
    #Assign merged dataframe back to the same name
    assign(df, merged_df)
    
  } else {
    warning(paste("Dataframe", df, "does not exist in the environment."))
  }
}




# ############################ #
# # How many isolates in modelling?# 
# ############################ #
df_names_modelling_test <- c("abaumanii_aminoglycosides_data", "abaumanii_carbapenems_data", 
              "abaumanii_fluoroquinolones_data", "escherichiacoli_carbapenems_data", 
              "escherichiacoli_cephalosporins_data", "klebsiellapneumoniae_carbapenems_data", 
              "klebsiellapneumoniae_cephalosporins_data", "pseudomonasaeroginosa_aminoglycosides_data", 
              "pseudomonasaeroginosa_cephalosporins_data")
# Initialize an empty list to store total sums
total_counts <- list()

# Loop through each dataframe and sum the "Total" field
for (df in df_names_modelling_test) {
  if (exists(df)) {
    temp_df <- get(df)
    # Sum the "Total" field
    total_sum <- sum(temp_df$Total, na.rm = TRUE)
    
    # Store the result in the list
    total_counts[[df]] <- total_sum
  } else {
    warning(paste("Dataframe", df, "does not exist in the environment."))
  }
}

# Convert to a dataframe for easier visualization
total_counts_df <- data.frame(
  Bug_Drug = names(total_counts),
  Total_Sum = unlist(total_counts)
)

# Print the result
print(total_counts_df)




# ############################ #
# # Scaling Data # 
# ############################ #



for (df in df_names) {
  if (exists(df)) {
    # Get the dataframe
    temp_df <- get(df)
    print(colnames(temp_df))
    # Check if the dataframe has at least 184 columns
    if (ncol(temp_df) >= 184) {
      # Scale columns 7 to 184
      temp_df[, 7:184] <- lapply(temp_df[, 7:184], scale)
      
      # Convert matrices back to numeric vectors
      temp_df[, 7:184] <- lapply(temp_df[, 7:184], function(x) as.numeric(x[, 1]))
      
      # Assign the transformed dataframe back to the original name
      assign(df, temp_df)
      
      # Print success message
      print(paste(df, "scaled and updated successfully."))
    } else {
      warning(paste("Skipping", df, "as it has less than 184 columns."))
    }
  } else {
    warning(paste("Dataframe", df, "does not exist in the environment."))
  }
}






# ############################ #
# # Splitting into Modelling and Prediction # 
# ############################ #



# Function to check and remove levels with only one observation
fix_factor_levels <- function(data, factor_vars) {
  for (var in factor_vars) {
    if (var %in% colnames(data)) {
      data[[var]] <- as.factor(data[[var]])  # Ensure it's a factor
      if (nlevels(data[[var]]) < 2) {
        warning(paste("Factor", var, "has <2 levels. Removing it from the model."))
        data <- data %>% select(-all_of(var))
      }
    }
  }
  return(data)
}

  
# Function to check and remove collinear variables using VIF
remove_collinear_features <- function(data, response_var) {
  # Fit a preliminary model to compute VIF
  lm_model <- lm(as.formula(paste(response_var, "~ .")), data = data)
  
  # Compute VIF values
  vif_values <- vif(lm_model)
  
  # Identify features with VIF > 5 (high collinearity)
  high_vif_features <- names(vif_values[vif_values > 5])
  
  # Remove high VIF features
  data <- data %>% select(-all_of(high_vif_features))
  
  return(data)
}

# Split the Data-frame into Modelling and Predictor Data Frames
for (df in df_names) {
  if (exists(df)) {
    # Get the dataframe
    temp_df <- get(df)
    
    # Ensure dataset has the required columns
    required_cols <- c("Resistant_per", "Total", "iso3", "Region")
    missing_cols <- setdiff(required_cols, colnames(temp_df))
    
    if (length(missing_cols) > 0) {
      warning(paste("Skipping", df, "- Missing columns:", paste(missing_cols, collapse = ", ")))
      next
    }
    
    # Apply standard filtering criteria to create initial `df_modelling`
    df_modelling <- temp_df[!is.na(temp_df$Resistant_per) & 
                              temp_df$Total >= 30 & 
                              !(temp_df$Resistant_per >= 0.99 | temp_df$Resistant_per <= 0.0005), ]
    
    # Apply GLM-based outlier detection **only on df_modelling**
    mad_results <- remove_outliers_mad(df_modelling, response_var = "Resistant_per",
                                       random_effects = c("iso3"))
    
    df_modelling <- mad_results$filtered_data  # Keep filtered data for modeling
    outlier_data <- mad_results$outliers       # Move outliers to df_prediction
    
    # Apply conditions for prediction dataset
    df_prediction <- temp_df[!( !is.na(temp_df$Resistant_per) & 
                                  temp_df$Total >= 30 & 
                                  !(temp_df$Resistant_per >= 0.99 | temp_df$Resistant_per <= 0.0005) ), ]
    
    # Add outliers to df_prediction
    df_prediction <- bind_rows(df_prediction, outlier_data)
    
    # Create unique names for modelling and prediction data
    modelling_name <- paste0(df, "_modelling")
    prediction_name <- paste0(df, "_prediction")
    
    # Assign them in the global environment
    assign(modelling_name, df_modelling, envir = .GlobalEnv)
    assign(prediction_name, df_prediction, envir = .GlobalEnv)
    
    # Print success message
    print(paste(df, "split into", modelling_name, "and", prediction_name, "successfully."))
  } else {
  
  }
}



# ############################ #
# # Splitting into Modelling into Test and Train# 
# ############################ #


# List of data-frame names
df_modelling_names <- c("abaumanii_aminoglycosides_data_modelling", "abaumanii_carbapenems_data_modelling", 
                        "abaumanii_fluoroquinolones_data_modelling", "escherichiacoli_carbapenems_data_modelling", 
                        "escherichiacoli_cephalosporins_data_modelling", "klebsiellapneumoniae_carbapenems_data_modelling", 
                        "klebsiellapneumoniae_cephalosporins_data_modelling", "pseudomonasaeroginosa_aminoglycosides_data_modelling", 
                        "pseudomonasaeroginosa_cephalosporins_data_modelling")


# Set seed for reproducibility
set.seed(123)

# Define output directory
train_test_output_dir <- "/Users/amr_user/Documents/Vivli Research Project/R Output/Testing_Training_datasets/"

# Create folder if it doesn't exist
if (!dir.exists(train_test_output_dir)) {
  dir.create(train_test_output_dir, recursive = TRUE)
}

# Initialize dataframe to store p-values
p_values_df <- data.frame(
  Bug_Drug = character(),
  P_Value_Resistant_per = numeric(),
  P_Value_Income_Group = numeric(),
  P_Value_Region = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each dataset and perform balanced train-test split
for (df in df_modelling_names) {
  if (exists(df)) {
    # Get the dataframe
    temp_df <- get(df)
    
    # Print dataset name and its columns for debugging
    print(paste("Processing:", df))
    print(colnames(temp_df))
    
    # Ensure dataset has the necessary columns
    required_columns <- c("Region", "Income group", "Resistant_per")
    missing_columns <- setdiff(required_columns, colnames(temp_df))
    
    if (length(missing_columns) > 0) {
      warning(paste("Skipping", df, "- Missing columns:", paste(missing_columns, collapse = ", ")))
      next
    }
    
    ### **1️⃣ Stratified Sampling (Ensuring Balance in "Region" & "Income group")** ###
    set.seed(123)  # Ensure reproducibility
    
    train_indices <- createDataPartition(temp_df$`Income group`, p = 0.75, list = FALSE)
    df_train <- temp_df[train_indices, ]
    df_test <- temp_df[-train_indices, ]
    
    ### **2️⃣ Ensure All "Income group" & "Region" Levels Exist in Both Train & Test** ###
    # Convert categorical variables to factors
    df_train$`Income group` <- as.factor(df_train$`Income group`)
    df_test$`Income group` <- as.factor(df_test$`Income group`)
    df_train$Region <- as.factor(df_train$Region)
    df_test$Region <- as.factor(df_test$Region)
    
    # Identify missing `Income group` levels in train and add them
    missing_income_groups <- setdiff(levels(df_test$`Income group`), levels(df_train$`Income group`))
    if (length(missing_income_groups) > 0) {
      print(paste("Adding missing Income groups to train dataset for:", df))
      extra_samples <- df_test %>% filter(`Income group` %in% missing_income_groups) %>% sample_n(1, replace = TRUE)
      df_train <- bind_rows(df_train, extra_samples)
    }
    
    # Identify missing `Region` levels in train and add them
    missing_regions <- setdiff(levels(df_test$Region), levels(df_train$Region))
    if (length(missing_regions) > 0) {
      print(paste("Adding missing Regions to train dataset for:", df))
      extra_samples <- df_test %>% filter(Region %in% missing_regions) %>% sample_n(1, replace = TRUE)
      df_train <- bind_rows(df_train, extra_samples)
    }
    
    # Ensure we have enough rows in both train and test sets
    if (nrow(df_train) < 2 | nrow(df_test) < 2) {
      warning(paste("Skipping", df, "- Not enough rows in train or test set. Train size:", nrow(df_train), "Test size:", nrow(df_test)))
      next
    }
    
    ### **3️⃣ Statistical Comparison of Train vs. Test** ###
    # Perform t-tests for continuous variables
    p_value_resistant_per <- tryCatch(t.test(df_train$Resistant_per, df_test$Resistant_per)$p.value, error = function(e) NA)
    
    # Convert "Income group" to numeric for statistical testing
    income_mapping <- c("Low income" = 1, "Lower middle income" = 2, "Upper middle income" = 3, "High income" = 4)
    df_train$Income_Numeric <- as.numeric(factor(df_train$`Income group`, levels = names(income_mapping), labels = income_mapping))
    df_test$Income_Numeric <- as.numeric(factor(df_test$`Income group`, levels = names(income_mapping), labels = income_mapping))
    
    # Perform Chi-square test for "Income group"
    income_train_table <- table(df_train$Income_Numeric)
    income_test_table <- table(df_test$Income_Numeric)
    p_value_income_group <- tryCatch(chisq.test(rbind(income_train_table, income_test_table))$p.value, error = function(e) NA)
    
    # Perform Chi-square test for "Region"
    region_train_table <- table(df_train$Region)
    region_test_table <- table(df_test$Region)
    p_value_region <- tryCatch(chisq.test(rbind(region_train_table, region_test_table))$p.value, error = function(e) NA)
    
    # Append results to dataframe
    new_row <- data.frame(
      Bug_Drug = df,
      P_Value_Resistant_per = p_value_resistant_per,
      P_Value_Income_Group = p_value_income_group,
      P_Value_Region = p_value_region,
      stringsAsFactors = FALSE
    )
    
    p_values_df <- rbind(p_values_df, new_row)
    
    # Print p-values for debugging
    print(paste(df, "- P-values:"))
    print(new_row)
    
    # Drop temporary Income_Numeric column
    df_train <- df_train %>% select(-Income_Numeric)
    df_test <- df_test %>% select(-Income_Numeric)
    
    # Create unique names for training and testing data
    train_name <- paste0(df, "_train")
    test_name <- paste0(df, "_test")
    
    # Assign the new datasets globally
    assign(train_name, df_train, envir = .GlobalEnv)
    assign(test_name, df_test, envir = .GlobalEnv)
    
    # Define file paths for saving
    train_file_path <- file.path(train_test_output_dir, paste0(train_name, ".xlsx"))
    test_file_path <- file.path(train_test_output_dir, paste0(test_name, ".xlsx"))
    
    # Save the train and test datasets as Excel files
    write_xlsx(df_train, train_file_path)
    write_xlsx(df_test, test_file_path)
    
    # Print success message
    print(paste(df, "split into", train_name, "and", test_name, "successfully."))
  } else {
    warning(paste("Dataframe", df, "does not exist in the environment."))
  }
}

# Save p-values dataframe
p_values_file_path <- file.path(train_test_output_dir, "Train_Test_P_Value_Comparison.xlsx")
write_xlsx(p_values_df, p_values_file_path)

print("Balanced train-test split completed.")











print(abaumanii_aminoglycosides_data_modelling_test)






# ############################ #
# # Child Models are trained here # 
# ############################ #

#dev.off() 

# set up cross-validation parameters for machine learning models in the caret package
cvCtrl <- trainControl(method = "cv", number = 5)


# Modelling Starts here..

# List of data-frame names
df_train_names <- c("abaumanii_aminoglycosides_data_modelling_train", "abaumanii_carbapenems_data_modelling_train", 
                    "abaumanii_fluoroquinolones_data_modelling_train", "escherichiacoli_carbapenems_data_modelling_train", 
                    "escherichiacoli_cephalosporins_data_modelling_train", "klebsiellapneumoniae_carbapenems_data_modelling_train", 
                    "klebsiellapneumoniae_cephalosporins_data_modelling_train", "pseudomonasaeroginosa_aminoglycosides_data_modelling_train", 
                    "pseudomonasaeroginosa_cephalosporins_data_modelling_train")


df_test_names <- c("abaumanii_aminoglycosides_data_modelling_test", "abaumanii_carbapenems_data_modelling_test", 
                    "abaumanii_fluoroquinolones_data_modelling_test", "escherichiacoli_carbapenems_data_modelling_test", 
                    "escherichiacoli_cephalosporins_data_modelling_test", "klebsiellapneumoniae_carbapenems_data_modelling_test", 
                    "klebsiellapneumoniae_cephalosporins_data_modelling_test", "pseudomonasaeroginosa_aminoglycosides_data_modelling_test", 
                    "pseudomonasaeroginosa_cephalosporins_data_modelling_test")




# Construct formula
sanitized_col_names <- colnames(abaumanii_aminoglycosides_data_modelling_train)[7:184]
sanitized_col_names <- paste0("`", sanitized_col_names, "`")
formula_model <- as.formula(paste("Resistant_per ~ Longitude + Latitude + year  +", paste(sanitized_col_names, collapse = " + ")))
# Print to check the formula
print(formula_model)



# Model 1 - GBM

# Initialize an empty list to store best R² values
rsquared_list <- list()
oos_results_df <- data.frame()

# Main training and evaluation loop
for (df_name in df_train_names) {
  if (exists(df_name)) {
    temp_df <- get(df_name)
    
    if (nrow(temp_df) > 1) {
      # Train GBM model
      gbm_model <- train(
        formula_model, 
        data = temp_df, 
        method = "gbm", 
        trControl = cvCtrl, 
        verbose = FALSE
      )
      
      # Save model to global env
      model_var_name <- paste0("gbm_model_", df_name)
      assign(model_var_name, gbm_model, envir = .GlobalEnv)
      
      # Feature importance summary
      model_summary_df <- as.data.frame(summary(gbm_model, plot = FALSE))
      file_name <- paste0(df_name, "_GBM_Summary.xlsx")
      file_path <- file.path("/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/gbm/", file_name)
      write_xlsx(model_summary_df, file_path)
      
      # In-sample R²
      best_r2 <- max(gbm_model$results$Rsquared)
      rsquared_list[[df_name]] <- best_r2
      
      # ----------------------------
      # Out-of-sample evaluation
      # ----------------------------
      test_df_name <- gsub("_train$", "_test", df_name)
      
      if (exists(test_df_name)) {
        test_df <- get(test_df_name)
        
        # Predict
        preds <- predict(gbm_model, newdata = test_df)
        actuals <- test_df$Resistant_per
        
        # Compute OOS R²
        ss_res <- sum((actuals - preds)^2)
        ss_tot <- sum((actuals - mean(actuals))^2)
        oos_r2 <- 1 - (ss_res / ss_tot)
        
        # Compute OOS RMSE
        oos_rmse <- sqrt(mean((actuals - preds)^2))
        
        # Append to result DF
        oos_results_df <- rbind(oos_results_df, data.frame(
          Model = df_name,
          OOS_R2 = round(oos_r2, 4),
          OOS_RMSE = round(oos_rmse, 4)
        ))
        
        print(paste("GBM model trained for", df_name, 
                    "| OOS R²:", round(oos_r2, 4), 
                    "| OOS RMSE:", round(oos_rmse, 4)))
      } else {
        warning(paste("Test dataframe", test_df_name, "does not exist. Skipping OOS evaluation."))
      }
      
    } else {
      warning(paste("Skipping", df_name, "as it has too few rows for training."))
    }
    
  } else {
    warning(paste("Dataframe", df_name, "does not exist in the environment."))
  }
}

# Save in-sample R²
rsquared_df <- data.frame(Model = names(rsquared_list), Best_Rsquared = unlist(rsquared_list))
write_xlsx(rsquared_df, "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/gbm/0.Training_GBM_Values.xlsx")

# Save OOS results
write_xlsx(oos_results_df, "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/gbm/0.OutOfSample_Performance_GBM.xlsx")



# Model 2 - Random Forest
# Define output directory for saving RF models

rf_output_dir <- "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/rf/"

# Create folder if it doesn't exist
if (!dir.exists(rf_output_dir)) {
  dir.create(rf_output_dir, recursive = TRUE)
}

# Initialize a list to store best R² values for RF
rsquared_list_rf <- list()
# Initialize a list for storing OOS performance
oos_results_rf_df <- data.frame()

# Loop through each dataframe and train the RF model
for (df_name in df_train_names) {
  if (exists(df_name)) {
    # Get the training dataframe
    temp_df <- get(df_name)
    
    if (nrow(temp_df) > 1) {
      #### Train Random Forest Model ####
      rf_model <- train(
        formula_model, 
        data = temp_df, 
        method = "rf", 
        trControl = cvCtrl, 
        verbose = FALSE
      )
      
      # Save RF model globally
      rf_model_var_name <- paste0("rf_model_", df_name)
      assign(rf_model_var_name, rf_model, envir = .GlobalEnv)
      
      #### Save Feature Importance ####
      rf_var_importance <- varImp(rf_model)
      rf_summary_df <- data.frame(
        Variable = rownames(rf_var_importance$importance),
        Importance = rf_var_importance$importance[, 1]
      )
      rf_file_path <- file.path(rf_output_dir, paste0(df_name, "_RF_Summary.xlsx"))
      write_xlsx(rf_summary_df, rf_file_path)
      
      #### In-sample R² ####
      best_r2_rf <- max(rf_model$results$Rsquared)
      rsquared_list_rf[[df_name]] <- best_r2_rf
      
      #### Out-of-Sample Evaluation ####
      test_df_name <- gsub("_train$", "_test", df_name)
      if (exists(test_df_name)) {
        test_df <- get(test_df_name)
        
        # Predict on test set
        preds_rf <- predict(rf_model, newdata = test_df)
        actuals_rf <- test_df$Resistant_per
        
        # Compute OOS R²
        ss_res_rf <- sum((actuals_rf - preds_rf)^2)
        ss_tot_rf <- sum((actuals_rf - mean(actuals_rf))^2)
        oos_r2_rf <- 1 - (ss_res_rf / ss_tot_rf)
        
        # Compute RMSE
        oos_rmse_rf <- sqrt(mean((actuals_rf - preds_rf)^2))
        
        # Save to result dataframe
        oos_results_rf_df <- rbind(oos_results_rf_df, data.frame(
          Model = df_name,
          OOS_R2_RF = round(oos_r2_rf, 4),
          OOS_RMSE_RF = round(oos_rmse_rf, 4)
        ))
        
        print(paste("RF model evaluated on test set for", df_name, 
                    "| OOS R²:", round(oos_r2_rf, 4), 
                    "| RMSE:", round(oos_rmse_rf, 4)))
      } else {
        warning(paste("Test dataframe", test_df_name, "does not exist. Skipping OOS evaluation."))
      }
    } else {
      warning(paste("Skipping", df_name, "as it has too few rows for training."))
    }
  } else {
    warning(paste("Dataframe", df_name, "does not exist in the environment."))
  }
}

# Save in-sample R²
rsquared_df_rf <- data.frame(
  Model = names(rsquared_list_rf), 
  Best_Rsquared_RF = unlist(rsquared_list_rf)
)
write_xlsx(rsquared_df_rf, file.path(rf_output_dir, "0.Training_RF_Values.xlsx"))

# Save out-of-sample performance
write_xlsx(oos_results_rf_df, file.path(rf_output_dir, "0.OutOfSample_Performance_RF.xlsx"))

# Final message
print("Random Forest model training and OOS evaluation completed.")







# Model 3 - XGBoost 
# Define output directory for saving XGBoost models
# Load required libraries

# Define output directory for saving XGBoost models
xgb_output_dir <- "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/xgboost/"
if (!dir.exists(xgb_output_dir)) {
  dir.create(xgb_output_dir, recursive = TRUE)
}

# Initialize result storage
rsquared_list_xgb <- list()
oos_results_xgb_df <- data.frame()

# Loop through each dataframe and train the XGBoost model
for (df_name in df_train_names) {
  if (exists(df_name)) {
    # Get the training dataframe
    temp_df <- get(df_name)
    
    if (nrow(temp_df) > 1) {
      #### Train XGBoost Model ####
      xgb_model <- tryCatch({
        train(
          formula_model, 
          data = temp_df, 
          method = "xgbTree", 
          trControl = cvCtrl, 
          verbose = FALSE
        )
      }, error = function(e) {
        warning(paste("Error in training XGBoost model for", df_name, ":", e$message))
        return(NULL)
      })
      
      # Skip if model training failed
      if (is.null(xgb_model) || any(is.na(xgb_model$results$Rsquared))) {
        warning(paste("XGBoost model failed for", df_name, "- skipping"))
        next
      }
      
      # Save XGBoost model globally
      xgb_model_var_name <- paste0("xgb_model_", df_name)
      assign(xgb_model_var_name, xgb_model, envir = .GlobalEnv)
      
      #### Save Feature Importance ####
      xgb_var_importance <- varImp(xgb_model)
      xgb_summary_df <- data.frame(
        Variable = rownames(xgb_var_importance$importance),  
        Importance = xgb_var_importance$importance[, 1]     
      )
      xgb_file_path <- file.path(xgb_output_dir, paste0(df_name, "_XGB_Summary.xlsx"))
      write_xlsx(xgb_summary_df, xgb_file_path)
      
      #### In-sample R² ####
      best_r2_xgb <- max(xgb_model$results$Rsquared, na.rm = TRUE)
      rsquared_list_xgb[[df_name]] <- best_r2_xgb
      
      #### Out-of-Sample Evaluation ####
      test_df_name <- gsub("_train$", "_test", df_name)
      if (exists(test_df_name)) {
        test_df <- get(test_df_name)
        
        preds_xgb <- predict(xgb_model, newdata = test_df)
        actuals_xgb <- test_df$Resistant_per
        
        ss_res_xgb <- sum((actuals_xgb - preds_xgb)^2)
        ss_tot_xgb <- sum((actuals_xgb - mean(actuals_xgb))^2)
        oos_r2_xgb <- 1 - (ss_res_xgb / ss_tot_xgb)
        oos_rmse_xgb <- sqrt(mean((actuals_xgb - preds_xgb)^2))
        
        oos_results_xgb_df <- rbind(oos_results_xgb_df, data.frame(
          Model = df_name,
          OOS_R2_XGB = round(oos_r2_xgb, 4),
          OOS_RMSE_XGB = round(oos_rmse_xgb, 4)
        ))
        
        print(paste("✅ XGBoost OOS evaluation for", df_name, 
                    "| R²:", round(oos_r2_xgb, 4), 
                    "| RMSE:", round(oos_rmse_xgb, 4)))
      } else {
        warning(paste("Test dataframe", test_df_name, "does not exist. Skipping OOS evaluation."))
      }
      
    } else {
      warning(paste("Skipping", df_name, "as it has too few rows for training."))
    }
  } else {
    warning(paste("Dataframe", df_name, "does not exist in the environment."))
  }
}

# Save in-sample R²
rsquared_df_xgb <- data.frame(
  Model = names(rsquared_list_xgb), 
  Best_Rsquared_XGB = unlist(rsquared_list_xgb)
)
write_xlsx(rsquared_df_xgb, file.path(xgb_output_dir, "0.Training_XGB_Values.xlsx"))

# Save out-of-sample performance
if (nrow(oos_results_xgb_df) > 0) {
  write_xlsx(oos_results_xgb_df, file.path(xgb_output_dir, "0.OutOfSample_Performance_XGB.xlsx"))
  print("📂 XGBoost out-of-sample performance saved successfully.")
} else {
 
}

print("🚀 XGBoost model training and OOS evaluation completed for all datasets.")





# Model 4 - GAM
# Load required libraries
library(mgcv)
library(writexl)

# Define output directory for saving GAM models
gam_output_dir <- "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/gam/"
if (!dir.exists(gam_output_dir)) {
  dir.create(gam_output_dir, recursive = TRUE)
}

# Initialize result containers
rsquared_list_gam <- list()
oos_results_gam_df <- data.frame()

# Loop through each dataframe and train the GAM model
for (df_name in df_train_names) {
  if (exists(df_name)) {
    temp_df <- get(df_name)
    
    # Sanitize column names to avoid errors
    colnames(temp_df) <- make.names(colnames(temp_df))
    
    # Construct GAM formula dynamically
    formula_gam_model <- as.formula(paste0(
      "Resistant_per ~ s(Longitude, Latitude) + ", 
      paste(
        colnames(temp_df)[!colnames(temp_df) %in% c("Resistant_per", "Longitude", "Latitude", "iso3", "Region", "Income.group", "Total")],
        collapse = " + "
      )
    ))
    
    print(formula_gam_model)
    
    #### Train GAM model ####
    gam_model <- gam(
      formula = formula_gam_model,
      data = temp_df,
      method = "REML"
    )
    
    # Save model to global environment
    gam_model_var_name <- paste0("gam_model_", df_name)
    assign(gam_model_var_name, gam_model, envir = .GlobalEnv)
    
    #### Extract Feature Importance ####
    parametric_terms <- summary(gam_model)$p.table[, "Estimate"]
    parametric_importance_df <- data.frame(
      Feature = names(parametric_terms),
      Importance = abs(parametric_terms)
    )
    
    smooth_terms <- summary(gam_model)$s.table[, "F"]
    smooth_importance_df <- data.frame(
      Feature = rownames(summary(gam_model)$s.table),
      Importance = smooth_terms
    )
    
    importance_df <- rbind(parametric_importance_df, smooth_importance_df)
    
    # Save importance
    gam_file_path <- file.path(gam_output_dir, paste0(df_name, "_GAM_Summary.xlsx"))
    write_xlsx(importance_df, gam_file_path)
    
    #### In-sample Pseudo-R² ####
    gam_r2 <- 1 - sum(residuals(gam_model)^2) / sum((temp_df$Resistant_per - mean(temp_df$Resistant_per))^2)
    rsquared_list_gam[[df_name]] <- gam_r2
    
    #### Out-of-sample Evaluation ####
    test_df_name <- gsub("_train$", "_test", df_name)
    if (exists(test_df_name)) {
      test_df <- get(test_df_name)
      colnames(test_df) <- make.names(colnames(test_df))  # sanitize again
      
      preds_gam <- predict(gam_model, newdata = test_df)
      actuals_gam <- test_df$Resistant_per
      
      ss_res_gam <- sum((actuals_gam - preds_gam)^2)
      ss_tot_gam <- sum((actuals_gam - mean(actuals_gam))^2)
      oos_r2_gam <- 1 - (ss_res_gam / ss_tot_gam)
      oos_rmse_gam <- sqrt(mean((actuals_gam - preds_gam)^2))
      
      oos_results_gam_df <- rbind(oos_results_gam_df, data.frame(
        Model = df_name,
        OOS_R2_GAM = round(oos_r2_gam, 4),
        OOS_RMSE_GAM = round(oos_rmse_gam, 4)
      ))
      
      print(paste("✅ GAM OOS evaluation for", df_name, 
                  "| R²:", round(oos_r2_gam, 4), 
                  "| RMSE:", round(oos_rmse_gam, 4)))
    } else {
      warning(paste("Test dataframe", test_df_name, "does not exist. Skipping OOS evaluation."))
    }
    
  } else {
    warning(paste("Dataframe", df_name, "does not exist in the environment."))
  }
}

# Save in-sample R²
rsquared_df_gam <- data.frame(
  Model = names(rsquared_list_gam), 
  Best_Rsquared_GAM = unlist(rsquared_list_gam)
)
write_xlsx(rsquared_df_gam, file.path(gam_output_dir, "0.Training_GAM_Values.xlsx"))

# Save out-of-sample performance
if (nrow(oos_results_gam_df) > 0) {
  write_xlsx(oos_results_gam_df, file.path(gam_output_dir, "0.OutOfSample_Performance_GAM.xlsx"))
  print("📂 GAM out-of-sample performance saved successfully.")
} else {

}

print("🚀 GAM model training and OOS evaluation completed for all datasets.")



# Load required libraries
library(caret)
library(writexl)

# Define output directory for saving SVM models
svm_output_dir <- "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/svm/"
if (!dir.exists(svm_output_dir)) {
  dir.create(svm_output_dir, recursive = TRUE)
}

# Initialize result containers
rsquared_list_svm <- list()
oos_results_svm_df <- data.frame()

# Loop through each training dataframe
for (df_name in df_train_names) {
  if (exists(df_name)) {
    temp_df <- get(df_name)
    
    if (nrow(temp_df) > 1) {
      #### Train SVM Model ####
      svm_model <- tryCatch({
        train(
          formula_model, 
          data = temp_df, 
          method = "svmRadial", 
          trControl = cvCtrl,
          verbose = FALSE
        )
      }, error = function(e) {
        warning(paste("Error in training SVM model for", df_name, ":", e$message))
        return(NULL)
      })
      
      if (is.null(svm_model) || any(is.na(svm_model$results$Rsquared))) {
        warning(paste("SVM model failed for", df_name, "- skipping"))
        next
      }
      
      # Save model to global environment
      svm_model_var_name <- paste0("svm_model_", df_name)
      assign(svm_model_var_name, svm_model, envir = .GlobalEnv)
      
      #### Feature Importance ####
      svm_var_importance <- varImp(svm_model)
      svm_summary_df <- data.frame(
        Variable = rownames(svm_var_importance$importance),
        Importance = svm_var_importance$importance[, 1]
      )
      svm_file_path <- file.path(svm_output_dir, paste0(df_name, "_SVM_Summary.xlsx"))
      write_xlsx(svm_summary_df, svm_file_path)
      
      #### In-sample R² ####
      best_r2_svm <- max(svm_model$results$Rsquared, na.rm = TRUE)
      rsquared_list_svm[[df_name]] <- best_r2_svm
      
      #### Out-of-Sample Evaluation ####
      test_df_name <- gsub("_train$", "_test", df_name)
      if (exists(test_df_name)) {
        test_df <- get(test_df_name)
        preds_svm <- predict(svm_model, newdata = test_df)
        actuals_svm <- test_df$Resistant_per
        
        ss_res_svm <- sum((actuals_svm - preds_svm)^2)
        ss_tot_svm <- sum((actuals_svm - mean(actuals_svm))^2)
        oos_r2_svm <- 1 - (ss_res_svm / ss_tot_svm)
        oos_rmse_svm <- sqrt(mean((actuals_svm - preds_svm)^2))
        
        oos_results_svm_df <- rbind(oos_results_svm_df, data.frame(
          Model = df_name,
          OOS_R2_SVM = round(oos_r2_svm, 4),
          OOS_RMSE_SVM = round(oos_rmse_svm, 4)
        ))
        
        print(paste("✅ SVM OOS evaluation for", df_name, 
                    "| R²:", round(oos_r2_svm, 4), 
                    "| RMSE:", round(oos_rmse_svm, 4)))
      } else {
        warning(paste("Test dataframe", test_df_name, "does not exist. Skipping OOS evaluation."))
      }
      
    } else {
      warning(paste("Skipping", df_name, "as it has too few rows for training."))
    }
  } else {
    warning(paste("Dataframe", df_name, "does not exist in the environment."))
  }
}

# Save in-sample R²
rsquared_df_svm <- data.frame(
  Model = names(rsquared_list_svm),
  Best_Rsquared_SVM = unlist(rsquared_list_svm)
)
write_xlsx(rsquared_df_svm, file.path(svm_output_dir, "0.Training_SVM_Values.xlsx"))

# Save out-of-sample performance
if (nrow(oos_results_svm_df) > 0) {
  write_xlsx(oos_results_svm_df, file.path(svm_output_dir, "0.OutOfSample_Performance_SVM.xlsx"))
  print("📂 SVM out-of-sample performance saved successfully.")
} else {
 
}

print("🚀 SVM model training and OOS evaluation completed for all datasets.")






# Model 6 - Cubist
# Define output directory for saving Cubist models
# Load necessary libraries
library(caret)
library(Cubist)
library(writexl)

# Define output directory
cubist_output_dir <- "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Child Models/cubist"
if (!dir.exists(cubist_output_dir)) {
  dir.create(cubist_output_dir, recursive = TRUE)
}

# Initialize result containers
rsquared_list_cubist <- list()
oos_results_cubist_df <- data.frame()

# Loop through each training dataset
for (df_name in df_train_names) {
  if (exists(df_name)) {
    temp_df <- get(df_name)
    
    if (nrow(temp_df) > 1) {
      # Sanitize column names
      colnames(temp_df) <- make.names(colnames(temp_df))
      
      # Create formula
      formula_cubist_model <- as.formula(paste0(
        "Resistant_per ~ ",
        paste(
          colnames(temp_df)[!colnames(temp_df) %in% c("Resistant_per", "iso3", "Income.group", "Income.group", "Region", "Total")],
          collapse = " + "
        )
      ))
      
      print(formula_cubist_model)
      
      # Train Cubist model
      cubist_model <- train(
        formula_cubist_model,
        data = temp_df,
        method = "cubist",
        trControl = cvCtrl
      )
      
      # Save model
      cubist_model_var_name <- paste0("cubist_model_", df_name)
      assign(cubist_model_var_name, cubist_model, envir = .GlobalEnv)
      
      # Variable importance
      cubist_var_importance <- varImp(cubist_model)
      cubist_summary_df <- data.frame(
        Variable = rownames(cubist_var_importance$importance),
        Importance = cubist_var_importance$importance[, 1]
      )
      
      # Save importance to Excel
      cubist_file_path <- file.path(cubist_output_dir, paste0(df_name, "_Cubist_Summary.xlsx"))
      write_xlsx(cubist_summary_df, cubist_file_path)
      
      # In-sample R²
      best_r2_cubist <- max(cubist_model$results$Rsquared, na.rm = TRUE)
      rsquared_list_cubist[[df_name]] <- best_r2_cubist
      
      #### Out-of-sample evaluation ####
      test_df_name <- gsub("_train$", "_test", df_name)
      if (exists(test_df_name)) {
        test_df <- get(test_df_name)
        colnames(test_df) <- make.names(colnames(test_df))
        
        preds_cubist <- predict(cubist_model, newdata = test_df)
        actuals_cubist <- test_df$Resistant_per
        
        ss_res_cubist <- sum((actuals_cubist - preds_cubist)^2)
        ss_tot_cubist <- sum((actuals_cubist - mean(actuals_cubist))^2)
        oos_r2_cubist <- 1 - (ss_res_cubist / ss_tot_cubist)
        oos_rmse_cubist <- sqrt(mean((actuals_cubist - preds_cubist)^2))
        
        oos_results_cubist_df <- rbind(oos_results_cubist_df, data.frame(
          Model = df_name,
          OOS_R2_Cubist = round(oos_r2_cubist, 4),
          OOS_RMSE_Cubist = round(oos_rmse_cubist, 4)
        ))
        
        print(paste("✅ Cubist OOS evaluation for", df_name, 
                    "| R²:", round(oos_r2_cubist, 4), 
                    "| RMSE:", round(oos_rmse_cubist, 4)))
      } else {
        warning(paste("Test dataframe", test_df_name, "does not exist. Skipping OOS evaluation."))
      }
      
    } else {
      warning(paste("Skipping", df_name, "as it has too few rows for training."))
    }
  } else {
    warning(paste("Dataframe", df_name, "does not exist in the environment."))
  }
}

# Save in-sample R²
rsquared_df_cubist <- data.frame(
  Model = names(rsquared_list_cubist),
  Best_Rsquared_Cubist = unlist(rsquared_list_cubist)
)
write_xlsx(rsquared_df_cubist, file.path(cubist_output_dir, "0.Training_Cubist_Values.xlsx"))

# Save out-of-sample performance
if (nrow(oos_results_cubist_df) > 0) {
  write_xlsx(oos_results_cubist_df, file.path(cubist_output_dir, "0.OutOfSample_Performance_Cubist.xlsx"))
  print("📂 Cubist out-of-sample performance saved successfully.")
} else {
 
}

print(" Cubist model training and OOS evaluation completed for all datasets.")




##Compute the best child models
library(readxl)
library(dplyr)
library(stringr)

# Define folder path
folder_path <- "/Users/amr_user/Documents/Vivli Research Project/R Output/ooskpis/"

# List all Excel files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Read, rename, and merge
combined_oos_df <- lapply(file_list, function(file) {
  df <- read_xlsx(file, col_names = TRUE)[, 1:3]
  colnames(df) <- c("Bug_Drug", "OOS_R2", "OOS_RMSE")  # Rename columns
  df$Model_Type <- str_extract(basename(file), "(?<=Performance_)[A-Za-z]+")  # Extract model name
  return(df)
}) %>%
  bind_rows()

# Preview
combined_oos_df <- combined_oos_df %>%
  mutate(OOS_R2 = ifelse(OOS_R2 < 0, 0.05, OOS_R2))


combined_oos_df <- combined_oos_df %>%
  mutate(
    Bug_Drug = str_replace(Bug_Drug, "_data_modelling_train$", "")
  )


print(head(combined_oos_df))

# Step 1: Compute R2 thresholds per Bug_Drug
r2_thresholds <- combined_oos_df %>%
  group_by(Bug_Drug) %>%
  summarise(
    R2_Max = max(OOS_R2, na.rm = TRUE),
    R2_SD  = sd(OOS_R2, na.rm = TRUE),
    R2_Threshold = R2_Max - R2_SD
  ) %>%
  ungroup()

# Step 2: Merge into combined_df and flag models to keep
combined_oos_df <- combined_oos_df %>%
  left_join(r2_thresholds, by = "Bug_Drug") %>%
  mutate(
    Keep_Model = ifelse(OOS_R2 >= R2_Threshold, TRUE, FALSE)
  )

# Preview the result
print(head(combined_oos_df))



selected_models_best_models_perbugdrug <- combined_oos_df %>%
  filter(Keep_Model == TRUE) %>%
  distinct(Bug_Drug, Model_Type)

# View the result
print(selected_models_best_models_perbugdrug)


# Optional: Save to file
# write.csv(combined_df, file = paste0(folder_path, "Combined_OOS_Metrics.csv"), row.names = FALSE)


# ############################ # # ############################ # # ############################ # 
# # Creating Stacked Child Models and training spatio-temporal Gaussian Process Regression # 
# ############################ # # ############################ # # ############################ #

# Loop through each dataframe and create stacked training datasets (only for selected models)
for (df_name in df_train_names) {
  if (exists(df_name)) {
    # Get training data
    train_data <- get(df_name)
    
    # Get base name (e.g., "abaumanii_aminoglycosides")
    base_name <- gsub("_data_modelling_train$", "", df_name)
    
    # Subset selected models for this bug-drug
    selected <- selected_models_best_models_perbugdrug %>% filter(Bug_Drug == base_name)
    
    # Initialize stacked dataframe with target and metadata
    stack_df_train <- data.frame(
      Resistant_per = train_data$Resistant_per,
      Year = train_data$year,
      Latitude = train_data$Latitude,
      Longitude = train_data$Longitude
    )
    
    # For GAM and Cubist: sanitize names
    train_data_sanitized <- train_data
    colnames(train_data_sanitized) <- make.names(colnames(train_data_sanitized))
    
    # Dynamically add model predictions only if in selected
    if ("GAM" %in% selected$Model_Type) {
      gam_model <- get(paste0("gam_model_", df_name))
      stack_df_train$gam <- predict(gam_model, newdata = train_data_sanitized, type = "response")
    }
    if ("GBM" %in% selected$Model_Type) {
      gbm_model <- get(paste0("gbm_model_", df_name))
      stack_df_train$gbm <- predict(gbm_model, newdata = train_data)
    }
    if ("RF" %in% selected$Model_Type) {
      rf_model <- get(paste0("rf_model_", df_name))
      stack_df_train$rf <- predict(rf_model, newdata = train_data)
    }
    if ("Cubist" %in% selected$Model_Type) {
      cubist_model <- get(paste0("cubist_model_", df_name))
      stack_df_train$cu <- predict(cubist_model, newdata = train_data_sanitized)
    }
    if ("XGBoost" %in% selected$Model_Type) {
      xgb_model <- get(paste0("xgb_model_", df_name))
      stack_df_train$xg <- predict(xgb_model, newdata = train_data)
    }
    if ("SVM" %in% selected$Model_Type) {
      svm_model <- get(paste0("svm_model_", df_name))
      stack_df_train$svm <- predict(svm_model, newdata = train_data)
    }
    
    # Save only if at least 1 model was used
    if (ncol(stack_df_train) > 5) {
      # Save stacked df globally
      stack_var_name <- paste0("stack_df_train_", base_name)
      assign(stack_var_name, stack_df_train, envir = .GlobalEnv)
      
      # Export to Excel
      write_xlsx(stack_df_train, paste0("/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Meta Models/stacked dataframes-training/stack_df_train_", base_name, ".xlsx"))
      
      print(paste(" Stacked training dataframe created for", base_name))
    } else {
      print(paste(" Skipping", base_name, "- no eligible models selected"))
    }
  }
}



# Loop through each dataframe and create stacked testing datasets
# Loop through each test dataframe and create stacked testing datasets
for (df_name in df_test_names) {
  if (exists(df_name)) {
    test_data <- get(df_name)
    
    # Extract the base name
    base_name <- gsub("_data_modelling_test$", "", df_name)
    
    # Get selected models for this bug-drug
    selected <- selected_models_best_models_perbugdrug %>% filter(Bug_Drug == base_name)
    
    # Initialize stacked dataframe with target and metadata
    stack_df_test <- data.frame(
      Resistant_per = test_data$Resistant_per,
      Year = test_data$year,
      Latitude = test_data$Latitude,
      Longitude = test_data$Longitude
    )
    
    # For GAM and Cubist: sanitize names
    test_data_sanitized <- test_data
    colnames(test_data_sanitized) <- make.names(colnames(test_data_sanitized))
    
    # Predict only for selected models
    if ("GAM" %in% selected$Model_Type) {
      gam_model <- get(paste0("gam_model_", base_name, "_data_modelling_train"))
      stack_df_test$gam <- predict(gam_model, newdata = test_data_sanitized, type = "response")
    }
    if ("GBM" %in% selected$Model_Type) {
      gbm_model <- get(paste0("gbm_model_", base_name, "_data_modelling_train"))
      stack_df_test$gbm <- predict(gbm_model, newdata = test_data)
    }
    if ("RF" %in% selected$Model_Type) {
      rf_model <- get(paste0("rf_model_", base_name, "_data_modelling_train"))
      stack_df_test$rf <- predict(rf_model, newdata = test_data)
    }
    if ("Cubist" %in% selected$Model_Type) {
      cubist_model <- get(paste0("cubist_model_", base_name, "_data_modelling_train"))
      stack_df_test$cu <- predict(cubist_model, newdata = test_data_sanitized)
    }
    if ("XGBoost" %in% selected$Model_Type) {
      xgb_model <- get(paste0("xgb_model_", base_name, "_data_modelling_train"))
      stack_df_test$xg <- predict(xgb_model, newdata = test_data)
    }
    if ("SVM" %in% selected$Model_Type) {
      svm_model <- get(paste0("svm_model_", base_name, "_data_modelling_train"))
      stack_df_test$svm <- predict(svm_model, newdata = test_data)
    }
    
    # Save only if model predictions are included
    if (ncol(stack_df_test) > 5) {
      stack_var_name <- paste0("stack_df_test_", base_name)
      assign(stack_var_name, stack_df_test, envir = .GlobalEnv)
      
      write_xlsx(
        stack_df_test,
        paste0("/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Meta Models/stacked dataframes-testing/stack_df_test_", base_name, ".xlsx")
      )
      
      print(paste(" Stacked test dataframe created for", base_name))
    } else {
      print(paste(" Skipping", base_name, "- no eligible models selected for test stacking"))
    }
  }
}



# ############################ # # ############################ # # ############################ # 
# # Training spatio-temporal Gaussian Process Regression # 
# ############################ # # ############################ # # ############################ #



stacked_training_names <- c("stack_df_train_abaumanii_aminoglycosides", "stack_df_train_abaumanii_carbapenems", 
              "stack_df_train_abaumanii_fluoroquinolones", "stack_df_train_escherichiacoli_carbapenems", 
              "stack_df_train_escherichiacoli_cephalosporins", "stack_df_train_klebsiellapneumoniae_carbapenems", 
              "stack_df_train_klebsiellapneumoniae_cephalosporins", "stack_df_train_pseudomonasaeroginosa_aminoglycosides", 
              "stack_df_train_pseudomonasaeroginosa_cephalosporins")



# Define output directories
meta_model_output_dir <- "/Users/amr_user/Documents/Vivli Research Project/R Output/Models/Meta Models/"
best_model_output_dir <- file.path(meta_model_output_dir, "Best GPR Models/")
dir.create(best_model_output_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize Lists to Store Performance Metrics
metrics_list <- list()
best_models_list <- list()
cv_results_list <- list()  # <--- ADD THIS LINE

# Loop through Each Stacked Training Dataset and Apply GPR
for (df_name in stacked_training_names) {
  print(df_name)
  
  # Extract the base name (e.g., "pseudomonasaeroginosa_cephalosporins")
  base_name <- gsub("^stack_df_train_", "", df_name)
  
  # Read Stacked Training Dataset
  stack_df_train <- read_excel(paste0(meta_model_output_dir, "stacked dataframes-training/stack_df_train_", base_name, ".xlsx"), sheet = "Sheet1")
  
  #### **Prepare Training Data** ####
  X_train <- stack_df_train %>%
    select(-Resistant_per) %>%
    mutate(across(everything(), as.numeric))  # Convert all columns to numeric
  
  y_train <- as.numeric(stack_df_train$Resistant_per)  # Convert target variable to numeric
  
  # Normalize Numeric Inputs
  normalize <- function(x) { (x - mean(x)) / sd(x) }
  X_train <- X_train %>% mutate(across(everything(), normalize))
  
  #### **2️⃣ Define Spatial-Temporal Groups** ####
  num_clusters <- 5  # Adjust based on dataset size
  spatial_clusters <- kmeans(X_train[, c("Latitude", "Longitude")], centers = num_clusters)$cluster
  X_train$Group_ID <- paste0(X_train$Year, "_", spatial_clusters)
  X_train$Group_ID <- as.factor(X_train$Group_ID)  # Convert to factor
  
  #### **3️⃣ Perform Group K-Fold Cross-Validation** ####
  k_folds <- 5
  cv_results <- list()
  best_rmse <- Inf  # Track best RMSE
  best_model <- NULL  # Store best model
  
  for (fold in 1:k_folds) {
    # Create Train-Test Split Using Group K-Fold
    train_indices <- which(X_train$Group_ID %in% levels(X_train$Group_ID)[-fold])
    test_indices <- which(X_train$Group_ID %in% levels(X_train$Group_ID)[fold])
    
    X_train_fold <- X_train[train_indices, ] %>% select(-Group_ID)
    y_train_fold <- y_train[train_indices]
    X_test_fold <- X_train[test_indices, ] %>% select(-Group_ID)
    y_test_fold <- y_train[test_indices]
    
    # Train GPR Model
    gpr_model <- km(
      formula = ~ .,  
      design = X_train_fold,  
      response = y_train_fold,  
      covtype = "matern5_2",
      control = list(parinit = "ARD")  
    )
    
    # Predict on Test Fold
    gpr_predictions <- predict(gpr_model, newdata = X_test_fold, type = "SK")
    
    # Compute Performance Metrics
    fold_rmse <- RMSE(as.vector(gpr_predictions$mean), y_test_fold)
    fold_r2 <- R2(as.vector(gpr_predictions$mean), y_test_fold)
    
    # Store Fold Results
    cv_results[[fold]] <- list(Fold = fold, RMSE = fold_rmse, R2 = fold_r2)
    
    # **Track Best Model Based on RMSE**
    if (fold_rmse < best_rmse) {
      best_rmse <- fold_rmse
      best_model <- gpr_model
    }
  }
  
  # Convert CV Results List to DataFrame
  cv_results_df <- do.call(rbind, lapply(cv_results, as.data.frame))
  
  # Save CV Results
  cv_results_list[[base_name]] <- cv_results_df
  
  # Store Best Performing Model Globally
  best_model_var_name <- paste0("best_gpr_model_", base_name)
  assign(best_model_var_name, best_model, envir = .GlobalEnv)
  best_models_list[[base_name]] <- list(Model_Name = best_model_var_name, Best_RMSE = best_rmse)
  
  # Save Best Model to File
  best_model_path <- paste0(best_model_output_dir, "best_gpr_model_", base_name, ".rds")
  saveRDS(best_model, best_model_path)
  
  print(paste("Best GPR model saved for", base_name, "with RMSE:", round(best_rmse, 4)))
  

}

# Convert CV Results to DataFrame and Save
final_cv_results_df <- do.call(rbind, lapply(names(cv_results_list), function(x) 
  cbind(Model = x, cv_results_list[[x]])))
write_xlsx(final_cv_results_df, paste0(meta_model_output_dir, "training-k-fold_performance.xlsx"))










# Define output directory for GPR results
gpr_results_output_dir <- file.path(meta_model_output_dir, "GPR_Test_Results/")
dir.create(gpr_results_output_dir, recursive = TRUE, showWarnings = FALSE)



# Set Number of Bootstrap Iterations
n_bootstrap <- 1000  # Adjust as needed

# Initialize a List to Store Bootstrapped Test Performance Metrics
bootstrapped_test_metrics <- list()

# Loop through Each Stacked Test Dataset and Apply Best GPR Model
for (df_name in df_test_names) {
  # Extract Base Name
  base_name <- gsub("_data_modelling_test$", "", df_name)
  print(base_name)
  
  # Read Stacked Testing Dataset
  stack_df_test <- read_excel(paste0(meta_model_output_dir, "stacked dataframes-testing/stack_df_test_", base_name, ".xlsx"), sheet = "Sheet1")
  
  # Retrieve the Corresponding Best Trained GPR Model
  gpr_model_var_name <- paste0("best_gpr_model_", base_name)
  
  #### **Prepare Test Data** ####
  X_test <- stack_df_test %>% select(-Resistant_per)  # Features
  y_test <- stack_df_test$Resistant_per               # Target variable
  
  # Normalize Numeric Inputs Using Training Statistics
  normalize <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }
  X_test <- X_test %>% mutate(across(everything(), normalize))
  
  # Load the Best Trained GPR Model
  gpr_model <- get(gpr_model_var_name)
  
  #### **Predict Using GPR**
  gpr_predictions <- predict(gpr_model, newdata = X_test, type = "SK")
  
  #### **Save GPR Test Predictions for Each Bug-Drug**
  gpr_test_results_df <- data.frame(
    Actual = y_test,  # Observed resistance values
    Predicted_Mean = gpr_predictions$mean,  # GPR Mean Predictions
    CI_Lower = gpr_predictions$mean - 1.96 * gpr_predictions$sd,  # 95% CI Lower Bound
    CI_Upper = gpr_predictions$mean + 1.96 * gpr_predictions$sd,  # 95% CI Upper Bound
    SD = gpr_predictions$sd,  # Standard Deviation for PDFs
    Bug_Drug = base_name,  # Bug-Drug Combination
    Type = "Test"  # Label Data Type
  )
  
  # Define file path
  gpr_test_results_path <- paste0(gpr_results_output_dir, "GPR_Test_Results_", base_name, ".xlsx")
  
  # Save to Excel
  write_xlsx(gpr_test_results_df, gpr_test_results_path)
  
  print(paste(" GPR test results saved for", base_name, "at", gpr_test_results_path))
  
  
  #### **Function for Bootstrapping Performance Metrics** ####
  boot_function <- function(data, indices) {
    # Sample data with replacement
    boot_sample <- data[indices, ]
    
    # Extract Bootstrapped X and Y
    X_boot <- boot_sample %>% select(-Actual_Resistant_per)
    y_boot <- boot_sample$Actual_Resistant_per
    
    # Apply GPR Model on Bootstrapped Sample
    boot_predictions <- predict(gpr_model, newdata = X_boot, type = "SK")
    
    # Compute Performance Metrics
    rmse_boot <- RMSE(as.vector(boot_predictions$mean), y_boot)
    r2_boot <- R2(as.vector(boot_predictions$mean), y_boot)
    mae_boot <- mean(abs(as.vector(boot_predictions$mean) - y_boot))
    pearson_corr_boot <- cor(as.vector(boot_predictions$mean), y_boot)
    nrmse_mean <- rmse_boot / mean(y_boot, na.rm = TRUE)
    # Return a named vector (ensures correct indexing)
    return(c(RMSE = rmse_boot, R2 = r2_boot, MAE = mae_boot, Pearson_Corr = pearson_corr_boot, nRMSE = nrmse_mean))
  }
  
  #### **Apply Bootstrapping** ####
  test_data_for_boot <- data.frame(Actual_Resistant_per = y_test, X_test)
  boot_results <- boot(data = test_data_for_boot, statistic = boot_function, R = n_bootstrap)
  
  # Convert Bootstrapped Results to DataFrame
  boot_results_df <- as.data.frame(boot_results$t)
  #colnames(boot_results_df) <- c("RMSE", "R2", "MAE", "Pearson_Corr")  # Ensure proper names
  colnames(boot_results_df) <- c("RMSE", "R2", "MAE", "Pearson_Corr", "nRMSE")
  # Extract Mean & 95% CI for Performance Metrics
  boot_means <- colMeans(boot_results_df, na.rm = TRUE)  # Mean values
  boot_ci <- apply(boot_results_df, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))  # 95% CI
  
  # Store Bootstrapped Metrics in a List
  # bootstrapped_test_metrics[[base_name]] <- c(
  #   RMSE_Mean = boot_means["RMSE"], RMSE_Lower95 = boot_ci["2.5%", "RMSE"], RMSE_Upper95 = boot_ci["97.5%", "RMSE"],
  #   R2_Mean = boot_means["R2"], R2_Lower95 = boot_ci["2.5%", "R2"], R2_Upper95 = boot_ci["97.5%", "R2"],
  #   MAE_Mean = boot_means["MAE"], MAE_Lower95 = boot_ci["2.5%", "MAE"], MAE_Upper95 = boot_ci["97.5%", "MAE"],
  #   Pearson_Corr_Mean = boot_means["Pearson_Corr"], Pearson_Corr_Lower95 = boot_ci["2.5%", "Pearson_Corr"], Pearson_Corr_Upper95 = boot_ci["97.5%", "Pearson_Corr"]
  # )
  bootstrapped_test_metrics[[base_name]] <- c(
    RMSE_Mean = boot_means["RMSE"], RMSE_Lower95 = boot_ci["2.5%", "RMSE"], RMSE_Upper95 = boot_ci["97.5%", "RMSE"],
    R2_Mean = boot_means["R2"], R2_Lower95 = boot_ci["2.5%", "R2"], R2_Upper95 = boot_ci["97.5%", "R2"],
    MAE_Mean = boot_means["MAE"], MAE_Lower95 = boot_ci["2.5%", "MAE"], MAE_Upper95 = boot_ci["97.5%", "MAE"],
    Pearson_Corr_Mean = boot_means["Pearson_Corr"], Pearson_Corr_Lower95 = boot_ci["2.5%", "Pearson_Corr"], Pearson_Corr_Upper95 = boot_ci["97.5%", "Pearson_Corr"],
    nRMSE_Mean = boot_means["nRMSE"], nRMSE_Lower95 = boot_ci["2.5%", "nRMSE"], nRMSE_Upper95 = boot_ci["97.5%", "nRMSE"]
  )
  
  # Print Bootstrapped Metrics
  print(paste("Bootstrapped Performance for", base_name, 
              "RMSE:", round(boot_means["RMSE"], 4), 
              "R²:", round(boot_means["R2"], 4),
              "MAE:", round(boot_means["MAE"], 4),
              "Pearson Correlation:", round(boot_means["Pearson_Corr"], 4)))
}

# Convert Bootstrapped Metrics List to DataFrame
bootstrapped_test_metrics_df <- do.call(rbind, lapply(names(bootstrapped_test_metrics), function(x) 
  cbind(Model = x, as.data.frame(t(bootstrapped_test_metrics[[x]])))))

# Save Bootstrapped Metrics to Excel
write_xlsx(bootstrapped_test_metrics_df, paste0(meta_model_output_dir, "bootstrapped_test_performance.xlsx"))

print(" Bootstrapping completed and results saved.")




# ############################ # # ############################ # # ############################ # 
# # Predicting Prevalence on unknown data now! # 
# ############################ # # ############################ # # ############################ #

print(abaumanii_aminoglycosides_data_prediction)

# List of data-frame names
df_prediction_names <- c("abaumanii_aminoglycosides_data_prediction", "abaumanii_carbapenems_data_prediction", 
                        "abaumanii_fluoroquinolones_data_prediction", "escherichiacoli_carbapenems_data_prediction", 
                        "escherichiacoli_cephalosporins_data_prediction", "klebsiellapneumoniae_carbapenems_data_prediction", 
                        "klebsiellapneumoniae_cephalosporins_data_prediction", "pseudomonasaeroginosa_aminoglycosides_data_prediction", 
                        "pseudomonasaeroginosa_cephalosporins_data_prediction")





# Ensure selected_models_df is in memory and has correct formatting
selected_models_best_models_perbugdrug <- selected_models_best_models_perbugdrug %>%
  mutate(Bug_Drug = gsub("_data_modelling_train", "", Bug_Drug))  # Clean base name

for (df_name in df_prediction_names) {
  
  base_name <- gsub("_data_prediction$", "", df_name)
  df_prediction <- get(df_name) %>% select(-Resistant_per, -Total)
  
  print(paste("🔄 Processing:", base_name))
  
  # Determine which models to keep for this bug-drug
  models_to_keep <- selected_models_best_models_perbugdrug %>%
    filter(Bug_Drug == base_name) %>%
    pull(Model_Type)

  # Prepare renamed prediction dataset for GAM/Cubist
  prediction_data_gam_cubist <- df_prediction
  colnames(prediction_data_gam_cubist) <- make.names(colnames(prediction_data_gam_cubist))

  # Initialize empty list to hold predictions
  pred_list <- list()

  # Predict using only selected models
  if ("GAM" %in% models_to_keep) {
    gam_model <- get(paste0("gam_model_", base_name, "_data_modelling_train"))
    pred_list$gam <- predict(gam_model, newdata = prediction_data_gam_cubist, type = "response")
  }
  if ("GBM" %in% models_to_keep) {
    gbm_model <- get(paste0("gbm_model_", base_name, "_data_modelling_train"))
    pred_list$gbm <- predict(gbm_model, newdata = df_prediction)
  }
  if ("RF" %in% models_to_keep) {
    rf_model <- get(paste0("rf_model_", base_name, "_data_modelling_train"))
    pred_list$rf <- predict(rf_model, newdata = df_prediction)
  }
  if ("Cubist" %in% models_to_keep) {
    cubist_model <- get(paste0("cubist_model_", base_name, "_data_modelling_train"))
    pred_list$cu <- predict(cubist_model, newdata = prediction_data_gam_cubist)
  }
  if ("XGBoost" %in% models_to_keep) {
    xg_model <- get(paste0("xgb_model_", base_name, "_data_modelling_train"))
    pred_list$xg <- predict(xg_model, newdata = df_prediction)
  }
  if ("SVM" %in% models_to_keep) {
    svm_model <- get(paste0("svm_model_", base_name, "_data_modelling_train"))
    pred_list$svm <- predict(svm_model, newdata = df_prediction)
  }

  # Add metadata
  pred_list$Year <- df_prediction$year
  pred_list$Latitude <- df_prediction$Latitude
  pred_list$Longitude <- df_prediction$Longitude

  # Combine into stacked dataframe
  stack_df_test <- as.data.frame(pred_list)

  print(paste("Successfully created stacked dataframe for", base_name))

  #### Predict using the best ST-GPR model ####
  gpr_model <- get(paste0("best_gpr_model_", base_name))

  # Normalize inputs
  normalize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  X_test <- stack_df_test %>% mutate(across(everything(), normalize))

  gpr_predictions <- predict(gpr_model, newdata = X_test, type = "SK")

  predictions_df <- data.frame(
    Resistant_per = pmax(0, pmin(1, as.vector(gpr_predictions$mean))), 
    CI_Lower = pmax(0, pmin(1, gpr_predictions$mean - 1.96 * gpr_predictions$sd)), 
    CI_Upper = pmax(0, pmin(1, gpr_predictions$mean + 1.96 * gpr_predictions$sd)),
    Year = df_prediction$year,
    Latitude = df_prediction$Latitude, 
    Longitude = df_prediction$Longitude,
    iso3 = df_prediction$iso3,
    `Income group` = df_prediction$`Income group`
  )

  # Write to file
  write_xlsx(predictions_df, paste0("/Users/amr_user/Documents/Vivli Research Project/R Output/Predictions/Predictions_", base_name, ".xlsx"))
}





# List of data-frame names
df_modelling_names <- c("abaumanii_aminoglycosides_data_modelling", "abaumanii_carbapenems_data_modelling", 
                         "abaumanii_fluoroquinolones_data_modelling", "escherichiacoli_carbapenems_data_modelling", 
                         "escherichiacoli_cephalosporins_data_modelling", "klebsiellapneumoniae_carbapenems_data_modelling", 
                         "klebsiellapneumoniae_cephalosporins_data_modelling", "pseudomonasaeroginosa_aminoglycosides_data_modelling", 
                         "pseudomonasaeroginosa_cephalosporins_data_modelling")

for (df_name in df_modelling_names) {
  
  # Extract base name by removing "_data_prediction"
  base_name <- gsub("_data_prediction$", "", df_name)
  #print(base_name)
  
  df_modelling <- get(df_name)
  
  print(base_name)
  
  ### **Store Predictions in a DataFrame** ####
  modelling_df <- data.frame(
    Resistant_per = df_modelling$Resistant_per, 
    Year = df_modelling$year,
    Latitude = df_modelling$Latitude, 
    Longitude = df_modelling$Longitude,
    iso3 = df_modelling$iso3,
    `Income group` = df_modelling$`Income group`
  )
  
  write_xlsx(modelling_df, paste0("/Users/amr_user/Documents/Vivli\ Research\ Project/R\ Output/Predictions/Knowndata_", base_name, ".xlsx") )
  
}








#Prediction Indicator Data
print(colnames(abaumanii_aminoglycosides_data))

# Load Required Libraries
library(tidyverse)
library(forecast)   # ARIMA
library(kernlab)    # GPR
library(mgcv)       # LOESS
library(readr)
library(writexl)

# Define Forecasting Years
forecast_years <- seq(2022, 2050)

# Extract Covariate Columns (7-184)
covariate_cols <- colnames(abaumanii_aminoglycosides_data)[7:184]

# Initialize a List to Store Forecasts
forecast_results <- list()

# Number of GPR Samples for Uncertainty Estimation
n_samples <- 100  # Can be increased for better uncertainty estimation

# Loop through Each Country (ISO3) and Region
for (iso3_code in unique(abaumanii_aminoglycosides_data$iso3)) {
  country_data <- abaumanii_aminoglycosides_data %>%
    filter(iso3 == iso3_code) %>%
    arrange(year) 
  
  region_name <- unique(country_data$Region)
  
  for (col in covariate_cols) {
    ts_data <- country_data %>% select(year, !!sym(col))
    
    # Check if Enough Data Exists for ARIMA
    if (sum(!is.na(ts_data[[col]])) > 10) {
      # Fit ARIMA Model
      arima_model <- auto.arima(ts_data[[col]], seasonal = FALSE)
      arima_forecast <- forecast(arima_model, h = length(forecast_years))
      
      # Extract Mean & SD from ARIMA
      forecast_arima_df <- data.frame(
        iso3 = iso3_code,
        Region = region_name,
        year = forecast_years,
        Variable = col,
        Forecast_ARIMA = as.numeric(arima_forecast$mean),
        SD_ARIMA = as.numeric((arima_forecast$upper[,2] - arima_forecast$lower[,2]) / 3.92)  # Approximate SD from 95% CI
      )
    } else {
      forecast_arima_df <- NULL
    }
    
    # Fit Gaussian Process Regression (GPR)
    gpr_model <- gausspr(as.matrix(ts_data$year), ts_data[[col]], kernel = "rbfdot")
    
    # Monte Carlo Sampling for GPR SD Estimation
    gpr_samples <- replicate(n_samples, {
      noise <- rnorm(length(forecast_years), mean = 0, sd = 0.1)  # Add small random noise
      predict(gpr_model, as.matrix(forecast_years)) + noise
    })
    
    # Compute Mean and Standard Deviation
    gpr_forecast_mean <- rowMeans(gpr_samples)
    gpr_forecast_sd <- apply(gpr_samples, 1, sd)  # Standard deviation of bootstrapped predictions
    
    # Store GPR Forecast
    forecast_gpr_df <- data.frame(
      iso3 = iso3_code,
      Region = region_name,
      year = forecast_years,
      Variable = col,
      Forecast_GPR = as.numeric(gpr_forecast_mean),
      SD_GPR = as.numeric(gpr_forecast_sd)  # Estimated Standard Deviation
    )
    
    # Apply LOESS Smoothing for Non-Stationary Covariates
    loess_model <- loess(ts_data[[col]] ~ ts_data$year, span = 0.3)
    loess_forecast <- predict(loess_model, newdata = forecast_years)
    
    # Store LOESS Forecast (No SD available for LOESS)
    forecast_loess_df <- data.frame(
      iso3 = iso3_code,
      Region = region_name,
      year = forecast_years,
      Variable = col,
      Forecast_LOESS = as.numeric(loess_forecast)
    )
    
    # Combine Forecasts
    combined_forecast <- full_join(forecast_arima_df, forecast_gpr_df, by = c("iso3", "Region", "year", "Variable")) %>%
      full_join(forecast_loess_df, by = c("iso3", "Region", "year", "Variable"))
    
    # Store Results
    forecast_results[[paste0(iso3_code, "_", col)]] <- combined_forecast
  }
}

# Convert Forecasts to a DataFrame
final_forecast_df <- bind_rows(forecast_results)

# Save Forecasts
write_csv(final_forecast_df, "/Users/amr_user/Documents/Vivli Research Project/R Output/Covariate_Forecasts_2022_2050.csv")
write_xlsx(final_forecast_df, "/Users/amr_user/Documents/Vivli Research Project/R Output/Covariate_Forecasts_2022_2050.xlsx")

# Print success message
print(" Forecasting completed and saved.")



