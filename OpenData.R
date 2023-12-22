library(data.table)
library(tidyr)
library(reshape2)
library(readxl)
library(purrr)
library(dplyr)
library(parallel)


library(countrycode)
add_continent_column <- function(dataframe, country_column_name = "Country") {
  
  country_vector <- as.character(dataframe[[country_column_name]])
  
  dataframe$Continent <- as.character(countrycode(country_vector, "country.name", "un.region.name"))
  
  cont_code <- countrycode(country_vector, "country.name", "un.region.code")
  #Because Latin American countries couldn't be done automatically it seems
  cont_code[country_vector %in% c("United States of America", "Canada", "Mexico")] <- 9
  
  dataframe$ContCode <- as.numeric(cont_code)
  
  return(dataframe)
}

read_and_process_excel <- function(file_path, countries) {
  data <- read_excel(file_path, col_types = "text")
  
  # Filter data for specified countries
  data <- filter(data, GeoAreaName %in% countries)
  
  return(data)
}

library(progress)

create_df_several_countries <- function(filepath, countries) {
  
  combined_data <- data.frame()
  total_countries <- length(countries)
  
  
  combined_data <- data.frame()  # Initialize an empty data frame

  data_list <- mclapply(filepath, read_and_process_excel,countries = countries)
  data_country <- bind_rows(data_list)
  combined_data <- bind_rows(combined_data, data_country)


  
  # Filter data for the specified country, location, sex, and age
  data <- combined_data %>%
    filter(is.na(`Location`) | `Location` == "ALLAREA") %>%
    filter(is.na(`Sex`) | `Sex` == "BOTHSEX") %>%
    filter(is.na(`Age`) | `Age` == "ALLAGE") %>%
    filter(is.na(`Type of speed`) | `Type of speed` == "ANYS") %>%
    filter(is.na(`Type of occupation`) | `Type of occupation` == "_T") %>%
    filter(is.na(`Migratory status`) | `Migratory status` == "_T") %>%
    filter(is.na(Activity) | Activity == "TOTAL" | Activity == "ISIC4_C") %>%
    filter(is.na(`Education level`) | `Education level` == "LOWSEC") %>%
    filter(Indicator != "12.4.1") %>%
    filter(is.na(`Type of waste treatment`) | `Type of waste treatment` == "INCINRT_EGY")
  # Select relevant columns
  data <- data %>%
    group_by(SeriesCode,GeoAreaName) %>%
    filter(var(Value) > 1e-6) %>%
    ungroup()
  data <-dplyr::select(data, TimePeriod, GeoAreaName, Indicator,SeriesCode,SeriesDescription , Value)
  data$Value <- as.numeric(data$Value)
  data$TimePeriod <- as.numeric(data$TimePeriod)
  data <- rename(data, Time = TimePeriod, Country = GeoAreaName)
  
  data <- data %>% filter(!is.na(Value))
  
  time_series_lengths <- data %>% group_by(SeriesCode) %>% summarize(series_length = n())
  
  #data <- data %>%
  #  inner_join(time_series_lengths, by = "SeriesCode") %>%
  #  filter(series_length >= 10) %>%
  #  dplyr::select(-series_length)
  
  data <- data %>%
    group_by(SeriesCode,Country) %>%
    filter(all(2010:2020 %in% Time)) %>%
    ungroup()%>%
    filter(Time %in% 2010:2020)
  
  
  series_counts <- data %>%
    group_by(SeriesCode) %>%
    summarize(NumCountries = n_distinct(Country))
  
  selected_series_codes <- series_counts %>%
    filter(NumCountries == total_countries) %>%
    pull(SeriesCode)
  
  data <- data %>%
    filter(SeriesCode %in% selected_series_codes)

  
  data$Name <- data$SeriesCode
  data$SeriesCode <- NULL
  

  
  data <-add_continent_column(data)
  data <- data %>%
    group_by(Name, Country) %>%
    mutate(Value_demeaned = Value - mean(Value, na.rm = TRUE)) %>%
    ungroup()
  data$Value <- data$Value_demeaned
  data <- data %>% dplyr::select(Name,Time,Value,Country,ContCode)
  

  return(data)
}
file_paths <- sprintf("C:\\Users\\kilia\\OneDrive\\Desktop\\Semesterthesis\\Data\\Goal%d.xlsx", 1:17)
specified_countries <- c("Switzerland","Japan","Germany","Republic of Korea","United States of America", "Canada", "Brazil")
#specified_countries <- c("Switzerland")
df <- create_df_several_countries(file_paths, specified_countries)
