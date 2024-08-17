setwd("~/vivli_amr")

library(dplyr)
library(lubridate)
library(lmtest)
library(sandwich)
library(zoo)
library(tseries)
library(tidyr)
library(ggplot2)
set.seed(123)

# Placeholder for final data storage
ccf_results <- list()

#antibiotic_data <- read.csv("data/blood_atlas.csv")
antibiotic_data <- read.csv("data/urine_atlas.csv")


#organism_antibiotic_pairs_blood <- list(
  #"Klebsiella pneumoniae" = c("Amikacin_I", "Amoxycillin.clavulanate_I", "Ampicillin_I", "Cefepime_I", "Ceftazidime_I", "Ceftriaxone_I", "Imipenem_I", "Levofloxacin_I", "Meropenem_I", "Minocycline_I", "Piperacillin.tazobactam_I", "Tigecycline_I", "Ampicillin.sulbactam_I", "Aztreonam_I", "Cefixime_I", "Ceftaroline_I", "Ceftazidime.avibactam_I", "Ciprofloxacin_I", "Colistin_I", "Doripenem_I", "Gentamicin_I", "Trimethoprim.sulfa_I", "Ceftolozane.tazobactam_I", "Meropenem.vaborbactam_I", "Cefpodoxime_I", "Ceftibuten_I"),
  #"Escherichia coli" = c("Amikacin_I", "Amoxycillin.clavulanate_I", "Ampicillin_I", "Cefepime_I", "Ceftazidime_I", "Ceftriaxone_I", "Imipenem_I", "Levofloxacin_I", "Meropenem_I", "Minocycline_I", "Piperacillin.tazobactam_I", "Tigecycline_I", "Ampicillin.sulbactam_I", "Aztreonam_I", "Cefixime_I", "Ceftaroline_I", "Ceftazidime.avibactam_I", "Ciprofloxacin_I", "Colistin_I", "Doripenem_I", "Gentamicin_I", "Trimethoprim.sulfa_I", "Ceftolozane.tazobactam_I", "Meropenem.vaborbactam_I", "Cefpodoxime_I", "Ceftibuten_I")
#)

organism_antibiotic_pairs <- list(
  "Klebsiella pneumoniae" = c("Amikacin_I", "Amoxycillin.clavulanate_I", "Ampicillin_I", "Cefepime_I", "Ceftazidime_I", "Ceftriaxone_I", "Imipenem_I", "Levofloxacin_I", "Meropenem_I", "Minocycline_I", "Piperacillin.tazobactam_I", "Tigecycline_I", "Ampicillin.sulbactam_I", "Aztreonam_I", "Cefixime_I", "Ceftaroline_I", "Ceftazidime.avibactam_I", "Ciprofloxacin_I", "Colistin_I", "Doripenem_I", "Gentamicin_I", "Trimethoprim.sulfa_I", "Ceftolozane.tazobactam_I", "Meropenem.vaborbactam_I", "Cefpodoxime_I", "Ceftibuten_I"),
  "Escherichia coli" = c("Amikacin_I", "Amoxycillin.clavulanate_I", "Ampicillin_I", "Cefepime_I", "Ceftazidime_I", "Ceftriaxone_I", "Imipenem_I", "Levofloxacin_I", "Meropenem_I", "Minocycline_I", "Piperacillin.tazobactam_I", "Tigecycline_I", "Ampicillin.sulbactam_I", "Aztreonam_I", "Cefixime_I", "Ceftaroline_I", "Ceftazidime.avibactam_I", "Ciprofloxacin_I", "Colistin_I", "Doripenem_I", "Gentamicin_I", "Trimethoprim.sulfa_I", "Ceftolozane.tazobactam_I", "Meropenem.vaborbactam_I", "Cefpodoxime_I", "Ceftibuten_I")
)

# Function to read and preprocess dataset 
read_and_filter_data <- function(organism, antibiotic_data) {
  df <- antibiotic_data %>% 
    mutate(gender = tools::toTitleCase(Gender),
           country = tools::toTitleCase(Country))
  # Filtering by organism
  df <- df %>% filter(Species == organism)
  df <- df %>% filter(Year > 2013)
  #df <- df %>% filter(Country == "Canada")
  
  return(df)
}

# Function to check if a series is stationary
check_stationary <- function(series) {
  series_p <- adf.test(series)$p.value
  return(series_p <= 0.05)
}

# Function to differ a series to make it stationary
my_diff <- function(x) {
  diff(x, differences = 1)
}

temp_ccf<- NULL
temp_ccf_pvalue <- NULL

make_stationary_and_calc_ccf <- function(series1, series2, anitibiotic_name1, anitibiotic_name2) {
  # print("Call to make_stationary")
  stat_series1 <- series1
  stat_series2 <- series2
  count<- 0
  
  if(length(stat_series1) != length(stat_series2)) {
    print("------ Diff Num rows detected -------------")
    #print(length(stat_series1))
    #print(length(stat_series2))
    sink()
    stop("Different number of rows detected")
  }
  
  while(TRUE) {
    s1_stationary <- check_stationary(stat_series1)
    s2_stationary <- check_stationary(stat_series2)
    
    if (count > 5 || is.na(s1_stationary) || is.na(s2_stationary))
      return(NaN)
    
    if (s1_stationary && s2_stationary)
      break
    
    stat_series1 <- my_diff(stat_series1)
    #print("Differenced TS1-----------------------------------")
    #print(stat_series1)
    stat_series2 <- my_diff(stat_series2)
    #print("Differenced TS2-----------------------------------")
    #print(stat_series2)
    count <- count + 1
    #print("Count is...........................")
    #print(count)

    }
  
  print(paste0(anitibiotic_name1,";" , anitibiotic_name2, ";", count))
  
  #print("From the function************************************************")
  #print(as.vector(stat_series1))
  #print(as.vector(stat_series2))
  
  ccf_values <- ccf(as.vector(stat_series1), 
                    as.vector(stat_series2), plot = "True", 
                    main = paste0(anitibiotic_name1, " ", anitibiotic_name2))
  ps <- 2 *(1 - pnorm(abs(ccf_values$acf), mean = 0, sd = 1/sqrt(ccf_values$n.used)))
  
  max_ccf_index = which.min(ccf_values$acf)
  max_ccf_p_value = ps[max_ccf_index]
  max_ccf_val = ccf_values$acf[max_ccf_index]
  
  temp_ccf <<- rbind(temp_ccf, c(anitibiotic_name1, 
                                 anitibiotic_name2,
                                 max_ccf_index-15,
                                 max_ccf_val,
                                 max_ccf_p_value))
}


for (organism in names(organism_antibiotic_pairs)) {
  antibiotics <- organism_antibiotic_pairs[[organism]]
  
  print(paste("Processing organism:", as.character(organism)))
  
  df <- read_and_filter_data(organism, antibiotic_data)
  
  # Create an empty data frame to store resistance percentages for this organism
  resistance_data <- data.frame()
  
  # Open the PDF file here for each organism
  pdf_filename <- paste0(as.character(organism), ".pdf")
  #pdf_filename <- "kleb_filtered_canada.pdf"
  pdf(pdf_filename, height = 6, width = 8)
  
  for (antibiotic in antibiotics) {
    if (!antibiotic %in% colnames(df)) {
      next
    }
    
    subset_data <- df[, c("Isolate.Id", "Year",  antibiotic)]
    
    print(paste("Processing antibiotic:", as.character(antibiotic)))
    
    # Group by Year and calculate resistance percentage
    
    ab_data <- subset_data %>%
      filter(.data[[antibiotic]] %in% c("Resistant", "Intermediate", "Susceptible")) %>%
      group_by(Year) %>%
      summarize(
        total_samples = n(),
        resistant_samples = sum(.data[[antibiotic]] == 'Resistant', na.rm = TRUE),
        resistance_percentage = (resistant_samples / total_samples) * 100
        
      )
    
    # Add antibiotic column name as a new column with resistance percentage values
    ab_data <- ab_data %>% 
      mutate(Antibiotic = antibiotic)
    
    # Reshape data so that years are rows and antibiotics are columns
    ab_data_wide <- ab_data %>%
      select(Year, Antibiotic, resistance_percentage) %>%
      pivot_wider(names_from = Antibiotic, values_from = resistance_percentage)
    
    # Combine with previous data
    resistance_data <- if (nrow(resistance_data) == 0) {
      ab_data_wide
    } else {
      full_join(resistance_data, ab_data_wide, by = "Year")
    }
  }
  
  # After processing all antibiotics for the organism
  if (!is.null(resistance_data) && ncol(resistance_data) > 1) {
    len_antibiotics <- length(antibiotics)
    for (i in 1:(len_antibiotics-1)) {
      for (j in (i + 1):len_antibiotics) {
        antibiotic1 <- antibiotics[i]
        antibiotic2 <- antibiotics[j]
        
        if (!(antibiotic1 %in% colnames(resistance_data)) || !(antibiotic2 %in% colnames(resistance_data))) {
          next
        }
        
        # Extract time series data for the two antibiotics
        series1 <- resistance_data[[antibiotic1]]
        #print("Original TS 1-----------------------------------")
        #print(series1)
        series2 <- resistance_data[[antibiotic2]]
        #print("Original TS 2-----------------------------------")
        #print(series2)
        
        # Check for NA values and skip if found
        if (any(is.na(series1)) || any(is.na(series2))) {
          next
        }
        
        #p <- ggplot(resistance_data, aes(x = Year)) +
          #geom_line(aes(y = series1, color = "Series1")) +
          #geom_line(aes(y = series2, color = "Series2")) +
          #labs(title = paste0(antibiotic1, " vs ", antibiotic2),
               #x = "Year",
               #y = "Resistance Percentage") +
         # scale_color_manual(values = c("Series1" = "blue", "Series2" = "red")) +
          #theme_minimal()
        #print (p)
        
        # Perform lead-lag analysis
        ccf_result <- make_stationary_and_calc_ccf(series1, series2, antibiotic1, antibiotic2)
        if (!is.na(ccf_result)) {
          ccf_results <- append(ccf_results, list(ccf_result))
          
        }
      }
    }
  }
  ccf_df <- do.call(rbind, lapply(ccf_results, as.data.frame))
  write.csv(ccf_df, paste0(as.character(organism), ".csv"), row.names = FALSE)
  # Close the PDF device after processing all pairs for the organism
  dev.off()
}



