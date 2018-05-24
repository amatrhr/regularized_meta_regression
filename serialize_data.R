# Generate a serialized list of data frames (RDS object) for shiny app ----
library(dplyr)
toy <- na.omit(read.csv("./toy_example.csv"))
toyY <- data.matrix(toy$y)
toyX <-  data.matrix(toy[,-1])
toyNames <- paste0(
  sample(LETTERS, nrow(toy), replace = TRUE),
  sample(month.name, nrow(toy), replace = TRUE))

murphy <- na.omit(read.csv("./Murphy.csv"))
murphyY <- data.matrix(murphy$Effect)
murphyX <- data.matrix(murphy[c(9:ncol(murphy))])
murphyNames <- murphy$Author

lightart <- na.omit(read.csv("./lightart.csv"))
lightartY <- data.matrix(lightart$Beta)
lightartX <- data.matrix(lightart[, 6:ncol(lightart)])
lightartNames <- lightart$Cohort

liu <- na.omit(read.csv("./liu_hair_morph_summ.csv"))
liu <- liu[liu$snp == "rs78924983",]
liuY <- data.matrix(liu$beta)
liuX <- data.matrix(liu[,c("purpose", "n", "mean_age", "sd_age", "min_age", "max_age",
  "pct_f", "pct_straight", "pct_wavy", "fea")])
liuNames <- data.matrix(liu$study)


ripke <- na.omit(read.csv("./ripke_et_al_st1_summ.csv"))
ripkeY <- data.matrix(ripke$lambda..pre.qc.)
ripkeX <- data.matrix(ripke %>% select(-starts_with("lambda"), -collection))
ripkeNames <- data.matrix(ripke$collection)

ripke_null <- na.omit(read.csv("./ripke_et_al_st2_summ.csv"))
ripke_nullY <- data.matrix(ripke_null$lambda..pre.qc.)
ripke_nullX <- data.matrix(ripke_null %>% select(-starts_with("lambda"), -collection))
ripke_nullNames <- data.matrix(ripke_null$collection)


shiny_data <- list(
  toy = list(X = toyX, Y = toyY, varnames = toyNames),
  murphy = list(X = murphyX, Y = murphyY, varnames = murphyNames),
  lightart = list(X = lightartX, Y = lightartY, varnames = lightartNames),
  liu = list(X = liuX, Y = liuY, varnames = liuNames),
  ripke = list(X = ripkeX, Y = ripkeY, varnames = ripkeNames),
  ripke_null = list(X = ripke_nullX, Y = ripke_nullY, varnames = ripke_nullNames)

  )

readr::write_rds(shiny_data, path = "./shiny_metareg_data.rds", compress = "gz")