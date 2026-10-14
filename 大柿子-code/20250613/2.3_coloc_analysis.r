# devtools::install_github("mrcieu/gwasglue")
# devtools::install_github("chr1swallace/coloc@main")
# remotes::install_github("stephenslab/susieR")
# devtools::install_github("boxiangliu/locuscomparer")
suppressMessages({
  library(dplyr)
  library(gwasglue)
  library(coloc)
  library(susieR)
  library(plinkbinr)
  library(doParallel)
  library(locuscomparer)
  library(RMySQL)
  library(glue)
  library(ggplot2)
  library(log4r)
}
)

pQTL_anno <- read.table("/work/database/MR/QTL/pQTL/WholeBlood_UKBPPP_2024/SMR/ID_Map_mito2.txt",sep = "\t",header = 1)
pQTL_id2map <- split(pQTL_anno$Gene,pQTL_anno$Probe) %>% lapply(function(x) unique(x))
eQTL_anno <- read.table("/work/database/MR/QTL/eQTL/eQTLGen_PhaseI/SMR/5e-8_1000/Id2map.txt",sep = "\t",header = 1,quote="")
eQTL_id2map <- split(eQTL_anno$Gene,eQTL_anno$Probe)
# mQTL_anno <- read.table("/work/database/MR/QTL/mQTL/WholeBlood_2024/EUR/SMR/5e-8_1000/Id2map.txt",sep = "\t",header = 1,quote="")
# mQTL_id2map <- split(mQTL_anno$Gene,mQTL_anno$Probe)

qtl_anno <- list(
  "eQTL" = eQTL_id2map,
  "pQTL" = pQTL_id2map
  # "mQTL" = mQTL_id2map
)

id2Trait <- read.csv("/work/run/songzx/MR/Project/MitoDys_SMR_CF/bin/ID2Trait.txt",sep = "\t",header = 1)
id2Trait_map <- split(id2Trait$FullName,id2Trait$ID)
## 函数区 ####
coloc_analysis <- function(dat_input,type,probe,gene,source){
  safe_runcoloc <- function(data_form1,data_form2,type,probe,gene,source) {  
    tryCatch({  
      my.res1 <- coloc.abf(dataset1 =data_form1 ,dataset2 = data_form2)
      dat_res <- data.frame(t(my.res1$summary)) %>% dplyr::mutate(probe = probe,
                                                                  Gene = gene,
                                                                  source = source,
                                                                  type = type )
    },  
    error = function(err) {  
      message("An error occurred: ", conditionMessage(err),"\n")  
      return(NULL)  
    })  
  }  
  
  dat_res <- NULL
  
  if(nrow(dat_input) > 10){
    dat_input <- dat_input %>%
      dplyr::group_by(SNP) %>%
      dplyr::slice_min(p_QTL,n=1,with_ties = FALSE) %>%
      dplyr::ungroup()
    dat_input[dat_input$p_QTL ==0,"p_QTL"] <- min(dat_input[dat_input$p_QTL != 0,"p_QTL"])
    dat_input[dat_input$p_GWAS ==0,"p_GWAS"] <-  min(dat_input[dat_input$p_GWAS != 0,"p_GWAS"])
    gwas1_form <- list(snp = dat_input$SNP,
                       position = dat_input$bp,
                       pvalues = dat_input$p_QTL,
                       # beta = dat_merge$beta_gwas1,
                       # varbeta = dat_merge$se_gwas1 ^ 2,
                       MAF = dat_input$MAF_QTL,
                       type = "quant",
                       N = dat_input$n_QTL)
    gwas2_form <-  list(snp = dat_input$SNP,
                        position =dat_input$bp,
                        beta = dat_input$b_GWAS,
                        varbeta = dat_input$se_GWAS^2,
                        type = "cc")
    
    dat_res <- safe_runcoloc(gwas1_form,gwas2_form,type,probe,gene,source)
  }
  return(dat_res)
}

susie_run <- function(dat_input,type,probe,gene,source){
  safe_runsusie <- function(data_form1,data_form2, maxit,repeat_until_convergence) {  
    tryCatch({  
      exp <- runsusie(data_form1,maxit = maxit,repeat_until_convergence =repeat_until_convergence,suffix = 1)
      out <- runsusie(data_form2,maxit = maxit,repeat_until_convergence =repeat_until_convergence,suffix = 2)
      susie.res <- coloc.susie(exp,out)
    },  
    error = function(err) {  
      message("An error occurred: ", conditionMessage(err))  
      return(NULL)  
    })  
  }  
  dat_res <- NULL

  if(nrow(dat_input) > 10){
    dat_input <- dat_input %>%
      group_by(SNP) %>%
      slice_min(p_QTL,n=1,with_ties = FALSE) %>%
      ungroup()
    LD <- ieugwasr::ld_matrix(dat_input$SNP, dat_input$SNP, with_alleles = F,bfile = "/work/database/MR/Reference/1000g/LDREF/EUR",
                              plink_bin =  plinkbinr::get_plink_exe())
    LD_snps <- dimnames(LD)[[1]]
    dat_input <- dat_input[dat_input$SNP %in% LD_snps, ]
    
    gwas1_form <- list(beta = dat_input$b_QTL, varbeta = dat_input$se_QTL^2,pvalues = dat_input$p_QTL,
                       LD = LD, snp = dat_input$SNP, type = "quant", MAF = dat_input$MAF_QTL, N = max(dat_input$n_QTL))
    
    gwas2_form <- list(snp = dat_input$SNP, beta = dat_input$b_GWAS, varbeta = dat_input$se_GWAS^2,
                       LD = LD, type = "cc", MAF = dat_input$MAF_GWAS, N = max(dat_input$n_GWAS))
    # gwas2_form <-  list(snp = dat_input$SNP,
    #                     position =dat_input$pos.outcome,
    #                     beta = dat_input$beta.outcome,
    #                     varbeta = dat_input$se.outcome^2,
    #                     type = "cc",
    #                     LD = as.matrix(dat_ld)[dat_input$SNP,dat_input$SNP],
    #                     N = max(dat_input$samplesize.outcome))
    susie.res <- safe_runsusie(gwas1_form,gwas2_form,maxit = 1000,repeat_until_convergence = FALSE)
    if(!is.null(susie.res$summary)){
      dat_res <- as.data.frame(susie.res$summary) %>% dplyr::mutate(probe = probe,
                                                                    Gene = gene,
                                                                    source = source,
                                                                    type = type )
    }
    return(dat_res)
  }else{
    return(NULL)
  }
}

parallel_func <- function(input_list,probe,gene,type,func,limit,packages){
  cl <- makeCluster(ifelse(length(input_list) > limit,limit,length(input_list)))
  registerDoParallel(cl)
  dat_para <- foreach(source = names(input_list),.combine = rbind,.packages = packages) %dopar% func(input_list[[source]],type,probe,gene,source)
  return(dat_para)
}
args <- commandArgs(T)
## 开始分析  ####
input_path <- "/work/run/songzx/MR/Project/MitoDys_SMR_CF/5_Coloc_res/"
file_id <- "CARDIoGRAM"
type <- "pQTL"
coloc_windows <- 100
suise_status <- 0

input_path <- args[1]
file_id <- args[2]
type <- args[3]
coloc_windows <- args[4]
suise_status <- args[5]
coloc_windows <- as.numeric(coloc_windows)

suise_status <- isTRUE(as.numeric(suise_status) == 1)


# 创建一个日志记录器  
logger <- create.logger(logfile =file.path(input_path,file_id,paste0(type,"_",coloc_windows,".log")), level = 'INFO')  

# 设置日志级别  
level(logger) <- 'INFO'  

input_path <- file.path(input_path,type,file_id,
                        sprintf("coloc_%s_%d.txt",type,coloc_windows))
dat_coloc <- data.table::fread(input_path,sep = "\t",data.table = F)
if(type != "mQTL"){
  dat_coloc$Gene <- sapply(dat_coloc$probe,function(x) qtl_anno[[type]][x]) %>% as.character()
}

## Coloc ####

dat_input_list <- lapply(plyr::dlply(dat_coloc,.variables = "Gene"),function(dat) 
  lapply(plyr::dlply(dat,.variables = "probe"),function(dat1) plyr::dlply(dat1,.variables = "source")))
info(logger,"正在进行初步共定位分析！")
dat_coloc_res <- lapply(names(dat_input_list),function(gene){
  lapply(names(dat_input_list[[gene]]),function(probe){
    parallel_func(dat_input_list[[gene]][[probe]],probe,gene,type,
                  coloc_analysis,30,packages = c("dplyr","coloc"))
  })
})
 
dat_coloc_res_all <- Reduce(rbind,Reduce(rbind,dat_coloc_res))
dat_coloc_res2 <- dat_coloc_res_all %>% na.omit()%>%
  dplyr::mutate(Type = ifelse(PP.H4.abf > 0.7,"Strong",
                              ifelse(PP.H4.abf > 0.5,"Mid","Weak")),
                method = "coloc")
info(logger,"初步共定位分析已完成！")
if(suise_status){
  info(logger,"正在进行Susie共定位分析！")
  dat_susie_res <- lapply(names(dat_input_list),function(probe) parallel_func(dat_input_list[[probe]],probe,source,type,
                                                                              susie_run,30, packages =c("dplyr","ieugwasr","plinkbinr","coloc","susieR")))
  dat_susie_res <- Reduce(rbind,dat_susie_res)
}else{
  dat_susie_res <- NULL
}

if(!is.null(dat_susie_res)){
  dat_susie_res2 <-  dat_susie_res %>%group_by(probe,source) %>%
    slice_max(PP.H4.abf,n=1,with_ties = FALSE) %>% ungroup() %>% dplyr::select(-c(hit1,hit2,idx1,idx2)) %>%
    mutate(Type = ifelse(PP.H4.abf > 0.8,"Strong",
                         ifelse(PP.H4.abf > 0.5,"Mid","Weak")),
           method = "susie")
  Result_All <- rbind(dat_coloc_res2,dat_susie_res2)
}else{
  Result_All <- dat_coloc_res2
}
info(logger,"Susie共定位分析已完成！")

Result_All_sig <- Result_All[Result_All$PP.H4.abf > 0.7,]
if(nrow(Result_All_sig) != 0){
  Result_All_sig$exposure <- paste0(Result_All_sig$Gene,"_",Result_All_sig$probe)
  Result_All_sig$outcome <- id2Trait_map[[file_id]]
}


write.table(Result_All,sprintf("%s/%s/Result_%s_%d_coloc.txt",input_path,file_id,type,coloc_windows),sep = "\t",quote = F,row.names = F)
write.table(Result_All_sig,sprintf("%s/%s/Result_%s_%d_coloc_sig.txt",input_path,file_id,type,coloc_windows),sep = "\t",quote = F,row.names = F)

info(logger,sprintf("共定位分析结果输出在%s!",sprintf("%s/%s/Result_%s_%d_coloc.txt",input_path,file_id,type,coloc_windows)))
## 共定位绘制点图 ####
# 如果两个结果都比较多的情况下才绘制这个图
# dat_coloc_plot <- coloc_use %>% 
#   dplyr::select('PP.H4.abf',exposure,method,outcome) %>%
#   tidyr::pivot_wider(names_from = 3,values_from = 1)
# 
# 
# ## locuscomparer绘图 ####
# protein_data_list <- list.files("6_coloc_data/pre/",pattern = glue::glue("Meta.*{best_windows}_coloc_pre.txt"),full.names = TRUE)
# coloc_res <- read.table("6_coloc_data/coloc_result_use.txt",sep = "\t",header = 1)
# 
# 
# for(rowid in 1:nrow(coloc_res)){
#   protein <- coloc_res[rowid,"exposure"]
#   trait <-coloc_res[rowid,"outcome"]
#   rsid <- coloc_res[rowid,"SNP"]
#   dat_plot_input <- protein_data_list[(grepl(paste0("_",protein,"_"),protein_data_list))&
#                                         (grepl(paste0("_",rsid,"_"),protein_data_list) &
#                                            (grepl(paste0("_",trait,"_"),protein_data_list)))]
#   dat_temp <- read.table(dat_plot_input,sep = "\t",header = 1)
#   in_fn1 = tempfile(fileext = ".ex1")
#   in_fn2 = tempfile(fileext = ".ex2")
#   write.table(dat_temp %>% dplyr::select(SNP,pval.outcome) %>% dplyr::rename("rsid"="SNP","pval"="pval.outcome") %>% distinct(),in_fn1,sep = "\t",quote = FALSE,row.names = FALSE)
#   write.table(dat_temp %>% dplyr::select(SNP,pval.exposure)%>% dplyr::rename("rsid"="SNP","pval"="pval.exposure") %>% distinct(),in_fn2,sep = "\t",quote = FALSE,row.names = FALSE)
#   plot <- locuscomparer::locuscompare(in_fn1 = in_fn1, in_fn2 = in_fn2, title1 = 'OA GWAS', title2 = paste0(protein,' pQTL'),
#                                       population = "EUR",combine = FALSE)
#   plot_total <- locuscomparer::locuscompare(in_fn1 = in_fn1, in_fn2 =in_fn2, title1 = 'OA GWAS', title2 = paste0(protein,' pQTL'),
#                                             population = "EUR",combine = TRUE)
#   ggsave(glue::glue("6_coloc_data/{trait}_{protein}_coloc_all.pdf"),plot = plot_total,width = 8,height = 6)
#   ggsave(glue::glue("6_coloc_data/{trait}_{protein}_coloc_two.pdf"),plot = plot[[1]],width = 8,height = 6)
#   ggsave(glue::glue("6_coloc_data/{trait}_{protein}_coloc_disease.pdf"),plot = plot[[2]],width = 8,height = 6)
#   ggsave(glue::glue("6_coloc_data/{trait}_{protein}_coloc_gene.pdf"),plot = plot[[3]],width = 8,height = 6)
#   message(glue("{protein} is over!"))
# }
# 
# write.csv(coloc_res,"result//Table3.csv",row.names = FALSE)
