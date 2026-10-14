rm(list = ls())
library(magrittr)

smr_res <- list.files("/work/run/songzx/MR/Project/MitoDys_SMR_DSZ/11_MultiSMR/Result",pattern = "*.smr",full.names = TRUE)
StepI_files <- smr_res[grepl("eQTL",smr_res)]
StepII_files <-  smr_res[grepl("pQTL",smr_res)]
dat_StepI <- Reduce(rbind,lapply(StepI_files,function(file){
  read.table(file,sep = "\t",header = TRUE)
})) %>% na.omit() %>%
  dplyr::select(1,3,17,18,19,20) %>%
  dplyr::mutate('OR(95% CI)' = paste0(
    sprintf("%.3f",exp(b_SMR)),"(",
    sprintf("%.3f",exp(b_SMR-1.96*se_SMR)),",",
    sprintf("%.3f",exp(b_SMR+1.96*se_SMR)),")"
  )) %>% dplyr::select(-c(3,4)) %>%
  dplyr::rename("p_StepI" = "p_SMR",
                "p_HEIDI_StepI" = "p_HEIDI",
                "OR(95% CI)_StepI" = "OR(95% CI)")

dat_StepII <- Reduce(rbind,lapply(StepII_files,function(file){
  read.table(file,sep = "\t",header = TRUE)
})) %>% na.omit() %>%
  dplyr::select(1,3,17,18,19,20) %>%
  dplyr::mutate('OR(95% CI)' = paste0(
    sprintf("%.3f",exp(b_SMR)),"(",
    sprintf("%.3f",exp(b_SMR-1.96*se_SMR)),",",
    sprintf("%.3f",exp(b_SMR+1.96*se_SMR)),")"
  )) %>% dplyr::select(-c(1,3,4)) %>%
  dplyr::rename("p_StepII" = "p_SMR",
              "p_HEIDI_StepII" = "p_HEIDI",
              "OR(95% CI)_StepII" = "OR(95% CI)")
dat_all <- dplyr::full_join(dat_StepI,dat_StepII) 
dat_sig_all <-   dplyr::full_join(dat_StepI %>%
                                    dplyr::filter(p_StepI < 0.05,
                                                  p_HEIDI_StepI > 0.01),
                                  dat_StepII %>% dplyr::filter(p_StepII < 0.05,
                                                               p_HEIDI_StepII > 0.01)) 
write.table(dat_all,"TwoStep_SMR.txt",sep = "\t",row.names = FALSE,quote = F)
write.table(dat_sig_all,"TwoStep_SMR_sig.txt",sep = "\t",row.names = FALSE,quote = F)
