library(tidyverse)

rna_pat_visits <- read.csv("~/OneDrive - Newcastle University/PhD/PPMI/rna_pat_visits.csv", header = TRUE)

bl <- read.csv("~/OneDrive - Newcastle University/PhD/PPMI/ppmi_bl_samples_filt.csv", header = T)

bl <- bl %>% 
  select(PATNO, status_clin, status_imag, clin_imag_agree, sex, age_bl, ppmi_enroll, filename_long)

rnaseq_availability <- merge(bl, rna_pat_visits, by = "PATNO")

write_csv(rnaseq_availability, "~/OneDrive - Newcastle University/PhD/PPMI/rnaseq_availability.csv")
