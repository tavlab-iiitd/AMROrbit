setwd("~/vivli_amr")

library(dplyr)
library(lubridate)
library(lmtest)
library(sandwich)
library(stats)
library(zoo)
set.seed(123)

#antibiotic_data <- read.csv("data/blood_atlas.csv")
antibiotic_data <- read.csv("data/urine_atlas.csv")

moving_avg <- function(x, n) {
  rollmean(x, k = n, fill = "extend", align = "center")
}

# Function to read and preprocess dataset 
read_and_filter_data <- function(organism, antibiotic_data) {
  df <- antibiotic_data %>% 
    mutate(gender = tools::toTitleCase(Gender),
           country = tools::toTitleCase(Country))
  # Filtering by organism
  df <- df %>% filter(Species == organism)
  df <- df %>% filter(Year > 2013)
  return(df)
}

#organism <- "Klebsiella pneumoniae"
organism <- "Escherichia coli"

antibiotics <- c ("Imipenem_I", "Meropenem_I", "Colistin_I")

for (antibiotic in antibiotics) {
  print(paste0("Processing......................", antibiotic))
  df <- read_and_filter_data(organism, antibiotic_data)
  subset_data <- df[, c("Isolate.Id", "Year", "Country", antibiotic)]
  
  accepted_vals <- c("Resistant", "Intermediate", "Susceptible")
  
  # Calculate the resistance percentage for each country, each year
  yearly_data <- subset_data %>%
    filter(.data[[antibiotic]] %in% accepted_vals) %>%
    group_by(Country, Year) %>%
    summarize(total_samples = n(),
              resistant_samples = sum(.data[[antibiotic]] == 'Resistant', na.rm = TRUE),
              resistance_percentage = (resistant_samples / total_samples) * 100)
  
  # Calculate global-level data
  global_data <- subset_data %>%
    filter(.data[[antibiotic]] %in% accepted_vals) %>%
    group_by(Year) %>%
    summarize(Country = "Global",  
              total_samples = n(),
              resistant_samples = sum(.data[[antibiotic]] == 'Resistant', na.rm = TRUE),
              resistance_percentage = (resistant_samples / total_samples) * 100)
  
  # Append global data to the existing yearly_data
  yearly_data_new <- rbind(yearly_data, global_data)
  
  yearly_data_new <- yearly_data_new %>% filter(!is.na(resistance_percentage))
  
  data_decomp <- yearly_data_new[, c("resistance_percentage", "Year", "Country")]
  
  # Identify start and end year and calculate the total number of data points
  start_year <- min(data_decomp$Year)
  end_year <- max(data_decomp$Year)
  total_years <- end_year - start_year + 1
  
  # Count the number of data points (years) for each country
  country_data_points <- data_decomp %>%
    group_by(Country) %>%
    summarize(DataPoints = n_distinct(Year))
  
  # Filter out countries with no missingness
  valid_countries <- country_data_points %>%
    filter(DataPoints == total_years) %>%
    pull(Country)
  
  data_decomp_final <- data_decomp %>%
    filter(Country %in% valid_countries)
  
  # Apply moving average to each group
  decomposed_data <- data_decomp_final %>%
    group_by(Country) %>%
    arrange(Year) %>%
    mutate(trend = moving_avg(resistance_percentage, n = 3)) %>%
    ungroup()
  
  # Filter out rows where trend is NA (because moving average will produce NAs at the boundaries)
  decomposed_data <- decomposed_data %>% filter(!is.na(trend))
  
  decomposed_data <- decomposed_data %>%
    mutate(seq_num = Year - start_year)
  
  decomposed_data$seq_num <- as.integer(decomposed_data$seq_num)
  decomposed_data$Country <- as.factor(decomposed_data$Country)
  decomposed_data$Country <- relevel(decomposed_data$Country, ref = "Global")
  
  model <- lm(trend ~ seq_num * factor(Country), data = decomposed_data)
  
  model_summary = summary(model)
  print(model_summary)
    
  robust_se <- coeftest(model, vcov = vcovHAC(model, type = "HC", cluster = "Country", group = decomposed_data$Country))
  
  print(robust_se)
  # Extract coefficients from the model summary
  coefficients <- model_summary$coefficients
  
  # Initialize empty vectors for slopes and intercepts
  intercepts <- numeric()
  slopes <- numeric()
  names <- character()
  
  # Separate global intercept and slope
  global_intercept <- coefficients["(Intercept)", "Estimate"]
  global_slope <- coefficients["seq_num", "Estimate"]
  
  # Add global values to the vectors
  names <- c(names, "Global")
  intercepts <- c(intercepts, global_intercept)
  slopes <- c(slopes, global_slope)
  
  # Extract country-specific intercepts and slopes
  for (name in rownames(coefficients)) {
    if (grepl("^factor\\(Country\\)", name)) {
      country_name <- sub("^factor\\(Country\\)", "", name)
      country_name <- sub("^:", "", country_name)  
      names <- c(names, country_name)
      intercept <- coefficients[name, "Estimate"]
      slope <- coefficients[paste0("seq_num:", name), "Estimate"]
      intercepts <- c(intercepts, intercept)
      slopes <- c(slopes, slope)
    }
  }
    
  # Create a data frame
  result_df <- data.frame(Country = names, intercept = intercepts, slope = slopes)
  
  # Save the result to a CSV file
  write.csv(result_df, paste0(paste0(organism, "_", antibiotic, "_urine_overall.csv")), row.names = FALSE)
  
}