# Load packages ----
library(tidyverse)
library(readxl)
library(stringr)

# read in data from files ----
st_1 <- read_xlsx("./Supplementary Table 1.xlsx", range = "B5:M14")
st_2 <- read_xlsx("./Supplementary Table 2.xlsx", range = "B5:AZ711")

# Clean names ----
st_1 <- rename(st_1, study = X__1, mean_age = Mean, sd_age = Sd, min_age = Min,
  max_age = Max, pct_F = `Female (%)`, pct_straight = `Straight (%)`,
  pct_wavy = `Wavy (%)`, pct_curly = `Curly (%)`) %>% select(-`_`)

st_1 <- mutate(st_1, Purpose = recode(Purpose, "D+P" = "D"),
  studypurpose = tolower(recode(study, "QIMR" = "disc", "RS" = "disc",
    "TwinsUK" = "disc")))


# First, a correspondence between study and appendix ----

studies <- 	c("ERF", "POL", "CANDELA", "US", "UYG", "TZL")

appendices <- c("__1", "__2", "__3", "__4", "__6", "__7")

rename_functions <- function(study, appendix){
  function(df){ df %>% rename_at( vars(ends_with(appendix)), funs(paste0(study,
    "_", gsub(., pattern = appendix, replacement = ""))))}
}

renamers <- map2(.x = studies, .y = appendices, .f = rename_functions)
names(renamers) <- paste0(studies, "_rename")

st_2 <- with(renamers, st_2 %>% ERF_rename() %>% TZL_rename() %>% UYG_rename()
  %>% CANDELA_rename() %>% POL_rename() %>% US_rename())

# Remove summary stats for non-asian samples and all samples
st_2 <- select(st_2, -contains("__"), -ERF_N, -N, -`_`)
st_2 <- rename(st_2, disc_fea = fEA, disc_beta = BETA, disc_se = SE, disc_p = P)

st_1 <- select_all(st_1, tolower)
st_2 <- select_all(st_2, tolower)

# Convert supplementary table 2 to long format ----
# 1. Gather all summary stats ---
# 2. Split into a column of study names and one of stat names ---
# 3. Spread ---


st_2 <- st_2 %>%
  gather(key = studystat, value = val, -loci, -snp, -chr, -bp,
  -ea, -nea) %>%
  separate(col = studystat, into = c("study", "stat")) %>%
  spread(stat, val)

# check the first snp, rs78924983, ----
st_2 %>% select(study, beta, fea, p, se, snp) %>% filter(snp == "rs78924983")
# check the last snp, rs366019, ----
st_2 %>% select(study, beta, fea, p, se, snp) %>% filter(snp == "rs366019")

# merge the two data frames ----
liu_hair_morph_summ <- left_join(st_1, st_2,
  by = c("studypurpose" = "study"))

# type conversion ----
liu_hair_morph_summ <- mutate(liu_hair_morph_summ,
  beta = as.numeric(beta),
  fea = as.numeric(fea),
  p = as.numeric(p),
  se = as.numeric(se)
  )


# initial analysis of this data (Sanity checks) ----
# check the first snp, rs78924983, ----
liu_hair_morph_summ %>% select(studypurpose, beta, fea, p, se, snp) %>%
  filter(snp == "rs78924983")
# check the last snp, rs366019, ----
liu_hair_morph_summ %>% select(studypurpose, beta, fea, p, se, snp, purpose) %>%
  filter(snp == "rs366019")


# write csv
write_csv(x = liu_hair_morph_summ, path = "./liu_hair_morph_summ.csv")
