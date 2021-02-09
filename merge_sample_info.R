library(tidyverse)

# set working directory
setwd("~/OneDrive - Newcastle University/PhD/PPMI")

# import existing sample file
samples <- read_csv("rnaseq_availability.csv", col_names = TRUE)

# import original cohort sample data to be added
og <- read.csv("PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020.csv", header = T)

# import prodromal cohort sample data to be added
prodromal <- read.csv("PPMI_Prodromal_Cohort_BL_to_Year_1_Dataset_Apr2020.csv", header = TRUE)

# merge og and prodromal datasets and write to file
merged_curated <- bind_rows(og, prodromal)
write_csv(merged_curated, file = "merged_curated.csv")

# merge datasets
merged_samples <- merge(samples, merged_curated, by = "PATNO")
write_csv(merged_samples, "merged_samples.csv")

# add genetic cohort and registry samples to all_samples
cohort_registry <- samples %>% 
  filter(status_clin == "GENPD" | status_clin == "GENUN" | status_clin == "REGPD" | status_clin == "REGUN")
write_csv(cohort_registry, "cohort_registry.csv")

## can't figure how to do it in r, so i'll do it in excel - add same variables from merged_samples to cohort_registry
cohort_registry_new <- read.csv("cohort_registry_new.csv", header = T)

# add REG and GEN rows
all_samples <- rbind(merged_samples, cohort_registry_new)
write_csv(all_samples, "all_samples.csv")


