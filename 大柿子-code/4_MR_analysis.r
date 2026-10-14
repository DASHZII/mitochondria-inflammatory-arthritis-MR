library(data.table)
library(TwoSampleMR)
library(dplyr)
library(stringr)
library(MRPRESSO)
library(forestploter)
library(cowplot)
library(ggplot2)
setwd("/work/run/songzx/MR/Project/MitoDys_SMR_CF/8_MR_Analysis/")
## 计算R2
calculate_r2 <- function(beta, eaf, samplesize, se) {
  numerator <- 2 * (beta^2) * eaf * (1 - eaf)
  denominator <- numerator + 2 * samplesize * eaf * (1 - eaf) * se^2
  return(numerator / denominator)
}
## 基础MR分析以及绘图
BaseMR <- function(MR_data,methods = c("mr_ivw","mr_simple_median","mr_weighted_median","mr_egger_regression")){
  heter_test <- 0.5
  presso <- NULL
  if(nrow(MR_data) == 1){
    methods <- c("mr_wald_ratio")
  }
  if(nrow(MR_data) > 3 ){
    presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = MR_data, NbDistribution = 1000,  
                        SignifThreshold = 0.05,seed = 123)
    if(presso$`MR-PRESSO results`$`Global Test`$Pvalue < 0.05){
      print(sprintf("有水平多效性，删去%s",paste(MR_data[presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,"SNP"],collapse = "/")))
      MR_data <- MR_data[-presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,]
      presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = MR_data, NbDistribution = 1000,  
                          SignifThreshold = 0.05,seed = 123)
    }
    
    heter_test <- presso$`MR-PRESSO results`$`Global Test`
  }
  
  
  
  res <- mr(MR_data,method_list = methods) %>% 
    mutate(OR = exp(b),CI_low =exp(b-1.96*se),CI_high = exp(b+1.96*se),
           method = ifelse(method== "Inverse variance weighted","IVW",method))
  
  heterogeneity <- mr_heterogeneity(MR_data)
  pleiotropy <- mr_pleiotropy_test(MR_data)
  
  if(!is.null(presso)){
    mr_presso <- cbind(res[1,c(1:4,6)],na.omit(presso$`Main MR results`) %>% mutate(OR = exp(`Causal Estimate`),
                                                                                    CI_low = OR-1.96*Sd,
                                                                                    CI_high = OR+1.96*Sd) %>%
                         dplyr::select(3:4,6:9) %>%
                         dplyr::rename("b"= "Causal Estimate",
                                       "se"="Sd",
                                       "pval" = "P-value") %>%
                         dplyr::mutate(method = "MR PRESSO"))
  }
  
  
  
  foreset_data <- res %>%
    select(3:ncol(.)) %>%
    mutate('OR(95% CI)' = paste0(format(OR,digits = 3),"(",format(CI_low,digits = 3),"-",format(CI_high,digits = 3),")")) %>%
    mutate(" "= paste(rep(" ", 12), collapse = " ")) 
  return(list(MR_data = MR_data,mr = res,heterogeneity = heterogeneity,presso = heter_test,pleiotropy=pleiotropy,foreset = foreset_data))
}


## MR分析绘图，包括散点图，漏斗图，留一法图
BasePlot <- function(res){
  res_single <- mr_singlesnp(res$MR_data)
  res_loo <- mr_leaveoneout(res$MR_data)
  
  p1 <- mr_scatter_plot(res$mr, res$MR_data)
  p2 <- mr_forest_plot(res_single)
  p4 <- mr_funnel_plot(res_single)
  p3 <- mr_leaveoneout_plot(res_loo)
  
  return(list(p1 = p1,p2 =p2,p3=p3,p4=p4))
}


input_path <- "/work/run/songzx/MR/Project/MitoDys_SMR_CF/raw_data"

dat_final_probe <- read.table("FInal_Sig_Probe.txt",sep = "\t",header = 1)

dat_meta <- read.table("/work/run/songzx/MR/Project/MitoDys_SMR_CF/MR_Analysis.txt",sep = "\t",header = 1)
disease_standard <-paste0(file.path(input_path,"Standard_GWAS",unique(dat_meta$File)),"_Standard_GWAS.txt")
names(disease_standard) <- c("Hypertension","Stroke")
cis_windows <- 500 ##http://dx.doi.org/10.1136/ard-2023-224459)
dat_gene_info <- data.table::fread("/work/run/songzx/Reference/Human/Genecode/HAVANA_gene_info.txt",sep = "\t")


id2Trait <- read.csv("/work/run/songzx/MR/Project/MitoDys_SMR_CF/bin/ID2Trait.txt",sep = "\t",header = 1)
id2Trait_map <- split(id2Trait$FullName,id2Trait$ID)

pQTL_files <- list.files("/work/run/songzx/MR/Project/MitoDys_SMR_CF/8_MR_Analysis/IVs/pQTL/",pattern = "*.txt",full.names = T)
eQTL_files <- list.files("/work/run/songzx/MR/Project/MitoDys_SMR_CF/8_MR_Analysis/IVs/eQTL/",pattern = "*.txt",full.names = T)
mQTL_files <- list.files("/work/run/songzx/MR/Project/MitoDys_SMR_CF/8_MR_Analysis/IVs/mQTL/",pattern = "*.txt",full.names = T)

# pQTL分析
dat_exp_pQTL <- Reduce(rbind,lapply(pQTL_files,function(inputFile){
  File_list <- stringr::str_split(basename(inputFile),"_",simplify = T)
  id_exp <- paste(File_list[1:(length(File_list)-2)],collapse = "_")
  exposure <- File_list[3]
  dat_exp <- data.table::fread(inputFile,data.table = F)
  dat_exp$id <- id_exp
  dat_exp$exposure <- exposure
  dat_exp_format <- format_data(dat_exp,type = "exposure",pos_col = "POS37",phenotype_col = "exposure")
})) %>%  dplyr::mutate(r2 = calculate_r2(beta.exposure,eaf.exposure,samplesize.exposure,se.exposure),
                       F = r2 * (samplesize.exposure - 2) / (1 - r2)) %>%
  dplyr::filter(F > 10)
dat_exp_pQTL <- dat_exp_pQTL %>% tidyr::separate("id.exposure",into = c("Probe","id.exposure"),sep ="_") %>%
  dplyr::left_join(dat_gene_info[,c("Probe","Chr","Pos","Gene")] %>%
                     dplyr::rename("exposure" = "Gene")) %>%
  dplyr::mutate(snp_type = ifelse(chr.exposure == Chr,ifelse(
    (pos.exposure > Pos -cis_windows*1000) &(pos.exposure < Pos +cis_windows*1000),"cis","trans"
  ),"trans"))


dat_out_pQTL <- Reduce(rbind,lapply(1:length(disease_standard),function(index){
  outcome_id <- gsub("_Standard_GWAS.txt","",basename(disease_standard[index]))
  dat_out <- data.table::fread(disease_standard[index],data.table = F)
  dat_out2 <- dat_out[dat_out$SNP %in% dat_exp_pQTL$SNP,]
  dat_out2$id <- outcome_id
  dat_out_format <-  format_data(dat_out2,type = "outcome",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "freq",
                                 beta_col = "b",pval_col = "p",samplesize_col = "n")
  dat_out_format$outcome <- names(disease_standard[index])
  return(dat_out_format)
}))

MR_data_pQTL <- harmonise_data(dat_exp_pQTL,
                          dat_out_pQTL,action=3)

MR_pQTL_Result <- Reduce(rbind,lapply(c("cis","trans"),function(type){
  dat_tmp <- plyr::dlply(MR_data_pQTL[MR_data_pQTL$snp_type == type,],c("id.exposure","id.outcome"),
              function(d){
                res <- NULL
                tryCatch({
                  res <- BaseMR(d)
                },error = function(e){
                  print(e)
                })
                plot <- NULL
                return(list(result = res,plot = plot))
              })
  
  Reduce(rbind,lapply(dat_tmp,function(d) d$result$mr)) %>%
    dplyr::mutate(mrType = type)
}))
pQTL_result <- plyr::dlply(MR_data_pQTL,c("id.exposure","id.outcome"),
                           function(d){
                             res <- NULL
                             tryCatch({
                               res <- BaseMR(d)
                             },error = function(e){
                               print(e)
                             })
                             plot <- NULL
                             return(list(result = res,plot = plot))
                           })

pQTL_mr_res <- Reduce(rbind,lapply(pQTL_result,function(d){
  d$result$mr
}))  %>%
  dplyr::mutate(mrType = "All")

MR_pQTL_Result <- rbind(MR_pQTL_Result,pQTL_mr_res)
pQTL_ivs_freqs <- table(MR_data_pQTL$exposure,MR_data_pQTL$snp_type) %>%
  as.data.frame() %>%
  tidyr::pivot_wider(id_cols = 1,values_from = 3,names_from = 2) %>%
  dplyr::mutate(All = cis+trans,
                QTL_type = "pQTL")

# eQTL分析  ####
dat_exp_eQTL <- Reduce(rbind,lapply(eQTL_files,function(inputFile){
  File_list <- stringr::str_split(basename(inputFile),"_",simplify = T)
  id_exp <- paste(File_list[1:(length(File_list)-2)],collapse = "_")
  exposure <-File_list[2]
  dat_exp <- data.table::fread(inputFile,data.table = F)
  dat_exp$id <- id_exp
  dat_exp$exposure <- exposure
  dat_exp_format <- format_data(dat_exp,type = "exposure",pos_col = "POS37",phenotype_col = "exposure")
})) %>%  dplyr::mutate(r2 = calculate_r2(beta.exposure,eaf.exposure,samplesize.exposure,se.exposure),
                       F = r2 * (samplesize.exposure - 2) / (1 - r2)) %>%
  dplyr::filter(F > 10)

dat_exp_eQTL <- dat_exp_eQTL %>% dplyr::mutate("Probe" = id.exposure) %>%
  dplyr::left_join(dat_gene_info[,c("Probe","Chr","Pos","Gene")] %>%
                     dplyr::rename("exposure" = "Gene")) %>%
  dplyr::mutate(snp_type = ifelse(chr.exposure == Chr,ifelse(
    (pos.exposure > Pos -cis_windows*1000) &(pos.exposure < Pos +cis_windows*1000),"cis","trans"
  ),"trans"))
eQTL_ivs_freqs <- table(dat_exp_eQTL$exposure,dat_exp_eQTL$snp_type) %>%
  as.data.frame() %>%
  tidyr::pivot_wider(id_cols = 1,values_from = 3,names_from = 2) %>%
  dplyr::mutate(All = cis+trans,
                QTL_type = "eQTL")

dat_out_eQTL <- Reduce(rbind,lapply(1:length(disease_standard),function(index){
  outcome_id <- gsub("_Standard_GWAS.txt","",basename(disease_standard[index]))
  dat_out <- data.table::fread(disease_standard[index],data.table = F)
  dat_out2 <- dat_out[dat_out$SNP %in% dat_exp_eQTL$SNP,]
  dat_out2$id <- outcome_id
  dat_out_format <-  format_data(dat_out2,type = "outcome",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "freq",
                                 beta_col = "b",pval_col = "p",samplesize_col = "n")
  dat_out_format$outcome <- names(disease_standard[index])
  return(dat_out_format)
}))

MR_data_eQTL <- harmonise_data(dat_exp_eQTL,
                               dat_out_eQTL,action=3)
MR_eQTL_Result <- Reduce(rbind,lapply(c("cis","trans"),function(type){
  dat_tmp <- plyr::dlply(MR_data_eQTL[MR_data_eQTL$snp_type == type,],c("id.exposure","id.outcome"),
                         function(d){
                           res <- NULL
                           tryCatch({
                             res <- BaseMR(d)
                           },error = function(e){
                             print(e)
                           })
                           plot <- NULL
                           return(list(result = res,plot = plot))
                         })
  
  Reduce(rbind,lapply(dat_tmp,function(d) d$result$mr)) %>%
    dplyr::mutate(mrType = type)
}))
eQTL_result <- plyr::dlply(MR_data_eQTL,c("id.exposure","id.outcome"),
                           function(d){
                             res <- NULL
                             tryCatch({
                               res <- BaseMR(d)
                             },error = function(e){
                               print(e)
                             })
                             plot <- NULL
                             return(list(result = res,plot = plot))
                           })

eQTL_mr_res <- Reduce(rbind,lapply(eQTL_result,function(d){
  d$result$mr
})) %>%
  dplyr::mutate(mrType = "All")

MR_eQTL_Result <- rbind(MR_eQTL_Result,eQTL_mr_res)
# mQTL分析  ####
dat_exp_mQTL <- Reduce(rbind,lapply(mQTL_files,function(inputFile){
  File_list <- stringr::str_split(basename(inputFile),"_",simplify = T)
  id_exp <- paste(File_list[1:(length(File_list)-2)],collapse = "_")
  exposure <- File_list[2]
  dat_exp <- data.table::fread(inputFile,data.table = F)
  dat_exp$id <- id_exp
  dat_exp$exposure <- exposure
  dat_exp_format <- format_data(dat_exp,type = "exposure",pos_col = "POS37",phenotype_col = "exposure")
})) %>%  dplyr::mutate(r2 = calculate_r2(beta.exposure,eaf.exposure,samplesize.exposure,se.exposure),
                       F = r2 * (samplesize.exposure - 2) / (1 - r2)) %>%
  dplyr::filter(F > 10)


# 使用检验过探针信息进行过滤
dat_exp_mQTL <- dat_exp_mQTL %>%  
  dplyr::left_join(dat_gene_info[,c("Probe","Chr","Pos","Gene")] %>%
                                     dplyr::rename("exposure" = "Gene")) %>%
  dplyr::mutate(snp_type = ifelse(chr.exposure == Chr,ifelse(
    (pos.exposure > Pos -cis_windows*1000) &(pos.exposure < Pos +cis_windows*1000),"cis","trans"
  ),"trans"))
dat_exp_mQTL <- dat_exp_mQTL %>% dplyr::filter((id.exposure %in% dat_final_probe$CpG) &(exposure %in% dat_final_probe$Exposure))

mQTL_ivs_freqs <-dat_exp_mQTL %>% dplyr::group_by(exposure,id.exposure,snp_type) %>%
  dplyr::summarise( n = n()) %>%
  tidyr::pivot_wider(names_from = 3,values_from = 4,values_fill = 0)%>%
  dplyr::mutate(All = cis+trans,
                QTL_type = "mQTL") %>%
  dplyr::rename("Var1" = "exposure",
                "CpG" = "id.exposure")

dat_out_mQTL <- dat_out_mQTL <- Reduce(rbind,lapply(1:length(disease_standard),function(index){
  outcome_id <- gsub("_Standard_GWAS.txt","",basename(disease_standard[index]))
  dat_out <- data.table::fread(disease_standard[index],data.table = F)
  dat_out2 <- dat_out[dat_out$SNP %in% dat_exp_mQTL$SNP,]
  dat_out2$id <- outcome_id
  dat_out_format <-  format_data(dat_out2,type = "outcome",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "freq",
                                 beta_col = "b",pval_col = "p",samplesize_col = "n")
  dat_out_format$outcome <- names(disease_standard[index])
  return(dat_out_format)
}))

MR_data_mQTL <- harmonise_data(dat_exp_mQTL,
                               dat_out_mQTL,action=3)

MR_mQTL_Result <- Reduce(rbind,lapply(c("cis","trans"),function(type){
  dat_tmp <- plyr::dlply(MR_data_mQTL[MR_data_mQTL$snp_type == type,],c("id.exposure","id.outcome"),
                         function(d){
                           res <- NULL
                           tryCatch({
                             res <- BaseMR(d)
                           },error = function(e){
                             print(e)
                           })
                           plot <- NULL
                           return(list(result = res,plot = plot))
                         })
  
  Reduce(rbind,lapply(dat_tmp,function(d) d$result$mr)) %>%
    dplyr::mutate(mrType = type)
}))

mQTL_result <- plyr::dlply(MR_data_mQTL,c("id.exposure","id.outcome"),
                           function(d){
                             res <- NULL
                             tryCatch({
                               res <- BaseMR(d)
                               },error = function(e){
                                 print(e)
                                      })
                             #plot <- BasePlot(res)
                             return(list(result = res))
                           })

mQTL_mr_res <- Reduce(rbind,lapply(mQTL_result,function(d){
  d$result$mr
})) %>% dplyr::mutate(mrType =  "All")
MR_mQTL_Result <-  rbind(MR_mQTL_Result,mQTL_mr_res)


MR_pQTL_Result$`OR(95% CI)` <- paste0(sprintf("%.2e",MR_pQTL_Result$OR),"(",sprintf("%.2e",MR_pQTL_Result$CI_low),",",sprintf("%.2e",MR_pQTL_Result$CI_high),")")
MR_eQTL_Result$`OR(95% CI)` <- paste0(sprintf("%.2e",MR_eQTL_Result$OR),"(",sprintf("%.2e",MR_eQTL_Result$CI_low),",",sprintf("%.2e",MR_eQTL_Result$CI_high),")")
MR_mQTL_Result$`OR(95% CI)` <- paste0(sprintf("%.2e",MR_mQTL_Result$OR),"(",sprintf("%.2e",MR_mQTL_Result$CI_low),",",sprintf("%.2e",MR_mQTL_Result$CI_high),")")

dat_pqtl_mr_res <- dat_meta %>% dplyr::rename("exposure" = "Gene",
                           "id.outcome"= "File") %>%
  dplyr::left_join(MR_pQTL_Result %>% dplyr::select(id.outcome,outcome,exposure,method,`OR(95% CI)`,pval,mrType))
dat_eqtl_mr_res <-  dat_meta %>% dplyr::rename("exposure" = "Gene",
                                               "id.outcome"= "File") %>%
  dplyr::left_join(MR_eQTL_Result %>% dplyr::select(id.outcome,outcome,exposure,method,`OR(95% CI)`,pval,mrType))
dat_mqtl_mr_res <-  dat_meta %>% dplyr::rename("exposure" = "Gene",
                                               "id.outcome"= "File") %>%
  dplyr::left_join(MR_mQTL_Result %>% dplyr::select(id.outcome,id.exposure,outcome,exposure,method,`OR(95% CI)`,pval,mrType))

ivs_Freq <- Reduce(rbind,list(mQTL_ivs_freqs,eQTL_ivs_freqs,pQTL_ivs_freqs))
ivs_Freq[is.na(ivs_Freq$CpG),"CpG"]  <- NA_character_
ivs_Freq <- dplyr::full_join(dat_meta %>% dplyr::select(1,2) ,
                 ivs_Freq %>%
                   dplyr::rename("Gene" = "Var1"))  %>%
  dplyr::mutate(Trait = id2Trait_map[File]  %>% as.character()) %>%
  dplyr::select(Trait,CpG,Gene,cis,trans,All,QTL_type)
write.table(ivs_Freq,"IVs_Freq_10000.txt",sep = "\t",row.names = F,quote = F)

rename_col <- c("OR(95% CI)","pval")

lapply(c("cis","trans","All"),function(type){
  dat_tmp <- dplyr::full_join(dat_mqtl_mr_res %>% 
                     dplyr::filter(method %in% c("IVW","Wald ratio")) %>%
                     dplyr::filter(mrType == type) %>%
                     dplyr::rename("OR(95% CI)_mQTL"="OR(95% CI)",
                                   "pval_mQTL"="pval",
                                   "CpG" = "id.exposure",
                                   "Exposure"="exposure",
                                   "File"="id.outcome",
                                   "Outcome" = "outcome") %>%
                     dplyr::right_join(dat_final_probe),
                   dat_eqtl_mr_res %>% 
                     dplyr::filter(method %in% c("IVW","Wald ratio"))%>% 
                     dplyr::filter(mrType == type) %>%
                     dplyr::rename("OR(95% CI)_eQTL"="OR(95% CI)",
                                   "pval_eQTL"="pval",
                                   "Exposure"="exposure",
                                   "File"="id.outcome",
                                   "Outcome" = "outcome") %>%
                     dplyr::right_join(dat_final_probe)) %>% 
    dplyr::full_join(
      dat_pqtl_mr_res %>% 
        dplyr::filter(method %in% c("IVW","Wald ratio"))%>% 
        dplyr::filter(mrType == type) %>%
        dplyr::rename("OR(95% CI)_pQTL"="OR(95% CI)",
                      "pval_pQTL"="pval",
                      "Exposure"="exposure",
                      "File"="id.outcome",
                      "Outcome" = "outcome") %>%
        dplyr::right_join(dat_final_probe)
    )
  write.table(dat_tmp,glue::glue("mr_analysis_res_{type}.txt"),sep = "\t",row.names = F,quote = F)
})

