rm(list = ls())
library(magrittr)
library(forestploter)


smr_res <- list.files("/work/run/songzx/MR/Project/MitoDys_SMR_DSZ/12_GTEx_tissue/SMR/Result/",pattern = "*.smr",full.names = TRUE)

dat_smr <- Reduce(rbind,lapply(smr_res,function(file){
  read.table(file,sep = "\t",header = TRUE) %>%
    dplyr::mutate(Tissue = strsplit(basename(file),"\\.")[[1]][1])
})) %>% na.omit() %>%
  dplyr::select(1,3,17,18,19,20,22) %>%
  dplyr::mutate('OR(95% CI)' = paste0(
    sprintf("%.3f",exp(b_SMR)),"(",
    sprintf("%.3f",exp(b_SMR-1.96*se_SMR)),",",
    sprintf("%.3f",exp(b_SMR+1.96*se_SMR)),")"
  ),
  or = sprintf("%.3f",exp(b_SMR)),
  ci_low = sprintf("%.3f",exp(b_SMR-1.96*se_SMR)),
  ci_hig = sprintf("%.3f",exp(b_SMR+1.96*se_SMR))) %>% dplyr::select(-c(3,4))
dat_smr$p_adj <- p.adjust(dat_smr$p_SMR,method = "bonferroni")

dat_smr$Gene <- "ATP5IF1"
table_s10 <- dat_smr[,c(5,2,6,3,4)]
write.table(table_s10,"TableS10.txt",sep = "\t",row.names = F)
figure_s10 <- dat_smr[,c(5,6,7,8,9,3)]
figure_s10$p_SMR <- sprintf("%.3e",figure_s10$p_SMR)
figure_s10 <- figure_s10 %>%
  dplyr::arrange(desc(or))
colnames(figure_s10)[6] <- "Pvalue"
tm <- forest_theme(base_size = 10,           # 设置基础字体大小
                   refline_col = "red4",     # 设置参考线颜色为红色
                   ci_col = "blue4",    # 设置箭头类型为闭合箭头
                   ci_pch = 15,
                   ci_Theight = 0.2)   # 设置脚注文字颜色为蓝色
p1 <- forest(figure_s10 %>%
               dplyr::select(1,2,6) %>%
               dplyr::mutate("  " = paste0(rep("    ",4),collapse = "")) %>%
               dplyr::select(1,2,4,3) ,
             est = figure_s10$or %>% as.numeric(),
             lower = figure_s10$ci_low%>% as.numeric(),
             upper = figure_s10$ci_hig%>% as.numeric(),
             ci_column = 3,
             ref_line = 1,
             xlim = c(0.5, 1.5),
             ticks_at = c(0.5, 1,1.5),
             theme = tm)

p1 <- add_border(p1,part = "header",row = 1)
p1
pdf("Figure8.pdf",width = 9,height = 6)
print(p1)
dev.off()
