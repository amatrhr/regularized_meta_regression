# Load libraries ----
library(readxl)
library(tidyverse)
library(stringr)
set.seed(56612)

# copy rename_functions from liu et al extraction,
# make separate  files later ---
rename_functions <- function(study, appendix){
  function(df){ df %>% rename_at( vars(ends_with(appendix)), funs(paste0(study,
  "_", gsub(., pattern = appendix, replacement = ""))))}
}
rename_Ripke_CC <- rename_functions(study = "controls", appendix = "__1")

# read in sheets 1, 2, 3, and 13 ----

st_1 <- read_xls("./ng.940-S2.xls", sheet = 1, range = "A3:S20", na = "n.a.")

st_2 <- read_xls("./ng.940-S2.xls", sheet = 2, range = "A3:S22", na = c("n.a.",
"n.a", "n..a."))

st_3_colnames <- names(read_xls("./ng.940-S2.xls", sheet = 3, range = "A2:S2"))
st_3 <- read_xls("./ng.940-S2.xls", sheet = 3, range = "A4:S20")

st_13_studies <- names(read_xls("./ng.940-S2.xls", sheet = 13, range = "B2:R2"))
st_13_methods <- read_xls("./ng.940-S2.xls", sheet = 13, range = "A3:A11", col_names = FALSE)
st_13_data <- read_xls("./ng.940-S2.xls", sheet = 13, range = "B3:R11", col_names = FALSE)

## unnest column names in sheet 1 ----
## Name, blank col, 10 case columns, blank col, 6 control columns ----

st_1 <- rename_Ripke_CC(st_1) %>% select(-controls_X, -contains("__"))

## unnest column names in sheet 2 ----
st_2 <- rename_Ripke_CC(st_2) %>% select(-controls_X, -contains("__"))

## assign column names to sheet 3 ----
colnames(st_3) <- st_3_colnames
st_3 <- select(st_3, -contains("__"))

## transpose sheet 13 and convert methods descriptors to short column names ----
st_13_data <- t(st_13_data)
first_last <- function(char_vec){
  map2(
  .x = lapply(char_vec,
    FUN = str_extract,
    pattern = regex("[:alpha:]*")),
  .y = lapply(char_vec,
    FUN = str_extract,
    pattern = regex("[:blank:][:alpha:]*$"))
  ,
  .f = function(x, y) {gsub(x = paste0(x,y), pattern = " ", replacement = "_")}
)
}
st_13_cols <- first_last(st_13_methods)

colnames(st_13_data) <- unlist(st_13_cols)

st_13 <- as.tibble(data.frame(Collection = st_13_studies, st_13_data))
st_13 <- mutate(st_13, Collection = as.character(Collection))


# parse sheet 1 ----
# clean years collected, ascertainment, instruments, exclusions ----
## Separate years_collected into: earliest_year, latest_year
## for cases and controls
## parse mean and sds of age at enrollment and age of onset
## into separate columns
parse_years_ripke <- function(dataframe) {

    dataframe %>%
      separate(`Years Collected`, c("first_cases", "last_cases"),
        convert = TRUE) %>%
      separate(`controls_Years Collected`, c("first_con", "last_con"),
        convert = TRUE) %>%
      separate(`Age at Onset`, c("mean_onset", "sd_onset"), convert = TRUE) %>%
      separate(`Age at Enrollment`, c("mean_enrol", "sd_enrol"),
        convert = TRUE) %>%
      separate(`controls_Age at Enrollment`, c("con_mean_enrol",
        "con_sd_enrol"), convert = TRUE)

}

st_1 <- parse_years_ripke(st_1)

# Case ascertainments, diagnoses, and diagnostic criteria ----
## recode ascertainment into inpatients, outpatients vs all others.
## recode diagnosis into SZ vs Other
## recode diagnostic criteria into DSM-IV vs other

st_1 <- st_1 %>%
    mutate(ascertainment = forcats::fct_lump(as.factor(Ascertainment)),
      diagnosis = forcats::fct_lump(as.factor(Diagnosis)),
    criteria = forcats::fct_lump(as.factor(`Diagnostic Criteria`))
)

st_1 <- st_1 %>% select(-Diagnosis, -Ascertainment, -`Diagnostic Criteria`)

# Instruments, exclusions, for cases/controls, ascertainment for controls ----
## Finally, to parse instruments and exclusions (And control ascertainment),
## extract unique strings from each variable
## make an indicator variable for each string

indicator_gen <- function(df, var, prefix = ""){
  ## List the unique, comma separated strings
  codes <- df[,var] %>% str_split(",") %>% unlist %>%
    str_replace_all(., pattern = "c\\(", replacement = "") %>%
    str_replace_all(., pattern ="[[:punct:]]" , replacement = "") %>%
    str_trim %>% unique()


  ## map: generate each indicator function
  ## apply indicator functions to the original variable
  ## store as a data frame to be bound to original df
  ## Shorten long strings
  ## don't include indicators which are true for only one observation

  indicator_df <- map_dfc(.x = codes,
    .f =  function(x){data.frame(grepl(x = unlist(data.frame(df[,var])),
      pattern =  x))}
    )

  codes[str_length(codes) > 8] <- first_last(codes[str_length(codes) > 8])
  colnames(indicator_df) <- paste0(prefix, codes)

  indicator_df <- indicator_df[,(colSums(indicator_df) > 1)]
  return(indicator_df)
}

st_1 <- bind_cols(st_1,
  indicator_gen(st_1, "Instruments", prefix = "i_"),
  indicator_gen(st_1, "Exclusions", prefix = "e_"),
  indicator_gen(st_1, "controls_Ascertainment", prefix = "con_a_"),
  indicator_gen(st_1, "controls_Instruments", prefix = "con_i_"),
  indicator_gen(st_1, "controls_Exclusions", prefix = "con_e_")
  ) %>% select(-Instruments, -Exclusions, -controls_Ascertainment,
    -controls_Instruments, -controls_Exclusions)

colnames(st_1) <- str_to_lower(colnames(st_1))

# parse sheet 2 ----
## Family history as double
st_2 <- mutate(st_2, `Family History` = as.double(`Family History`))
# years collected, ascertainment, instruments, exclusions ----
## Separate years_collected into: earliest_year, latest_year, duration

## Exclude UQ and ASRB from analysis
st_2 <- st_2 %>% filter((!grepl(x = Collection, pattern = "UQ")) )

## Convert Finnish studies' ages to NA
st_2 <- st_2 %>% mutate(
  `Age at Enrollment` = recode(`Age at Enrollment`,
    "born between 1940-1974" = NA_character_),
  `controls_Age at Enrollment` = recode(`controls_Age at Enrollment`,
    "20-70 years old" = NA_character_))

st_2 <- parse_years_ripke(st_2)

# Ascertainments, instruments, exclusions, for cases and controls ----

st_2 <- st_2 %>%  mutate(
  diagnosis = forcats::fct_lump(as.factor(Diagnosis)),
  criteria = forcats::fct_lump(as.factor(`Diagnostic Criteria`))
)

st_2 <- st_2 %>% select(-Diagnosis, -`Diagnostic Criteria`)


st_2 <- bind_cols(st_2,
  indicator_gen(st_2, "Ascertainment", prefix = "a_"),
  indicator_gen(st_2, "Instruments", prefix = "i_"),
  indicator_gen(st_2, "Exclusions", prefix = "e_"),
  indicator_gen(st_2, "controls_Ascertainment", prefix = "con_a_"),
  indicator_gen(st_2, "controls_Instruments", prefix = "con_i_"),
  indicator_gen(st_2, "controls_Exclusions", prefix = "con_e_")
) %>% select(-Ascertainment, -Instruments, -Exclusions, -controls_Ascertainment,
  -controls_Instruments, -controls_Exclusions)

colnames(st_2) <- str_to_lower(colnames(st_2))

# parse sheet 3 ----
colnames(st_3) <- str_to_lower(colnames(st_3))

# parse sheet 13 ----
colnames(st_13) <- str_to_lower(colnames(st_13))

st_13 <-  st_13 %>% mutate_at(vars(contains("_")),funs(recode(., "yes" = TRUE,
  "no" = FALSE)))

## Again, if only one study has a characteristic, drop the characteristic
morethanone <- function(column){
  return(as.logical((sum(column) > 1)*(sum(column) < length(column))))
}
st_13[,-1] <- st_13[,-1] %>% select_if(.predicate = morethanone)

# join sheets ----
sheet_1_3 <- left_join(st_1, st_3)
sheet_1_3_13 <- left_join(sheet_1_3, st_13)


sheet_2_3_13  <- bind_cols(st_2,
  st_3 %>% select(starts_with("lambda")) %>%
    sample_n(size = nrow(st_2), replace = TRUE),
  st_13 %>% select(contains("_")) %>%
    sample_n(size = nrow(st_2), replace = TRUE)
  )
# Sanity checks ----


# export csv file ----

## Export Stage 1 Studies' CSV with observed Lambda_GC
write_csv(sheet_1_3_13, path = "./ripke_et_al_st1_summ.csv")

## Export Stage 2 Studies' CSV with simulated Lambda_GCs and broken correlation
## between "Sheet 13" coding and indicator variables
write_csv(sheet_2_3_13, path = "./ripke_et_al_st2_summ.csv")