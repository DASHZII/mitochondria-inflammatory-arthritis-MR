library(magrittr)
args <- commandArgs(T)
##染色体长度处理
# dat_chrom_length <- read.table("/work/database/MR/Reference/chrom.size",sep = "\t")
# dat_chrom_length$V1 <- sapply(dat_chrom_length$V1,function(x) stringr::str_remove(x,"^chr"))
# dat_chrom_length <- dat_chrom_length[dat_chrom_length$V1 %in% seq(1,23),]
# chromLength_37 <- split(dat_chrom_length[dat_chrom_length$V3 == "GRCh37",]$V2,dat_chrom_length[dat_chrom_length$V3 == "GRCh37",]$V1)
# chromLength_38 <- split(dat_chrom_length[dat_chrom_length$V3 == "GRCh38",]$V2,dat_chrom_length[dat_chrom_length$V3 == "GRCh37",]$V1)
# chromLength <- list("GRCh37" = chromLength_37,
#      "GRCh38" = chromLength_38)
# qs::qsave(chromLength,"/work/database/MR/Reference/chrom.size.qs")
selCol <- c("SNP","chr","POS37","POS38","id","exposure")
# coloc_windows <- list(pQTL = 100,
#                       eQTL = 100,
#                       mQTL = 100)

mQTL_besd <- "/work/database/MR/QTL/mQTL/WholeBlood_2024/EUR/SMR/All"
eQTL_besd <- "/work/database/MR/QTL/eQTL/eQTLGen_PhaseI/SMR/All/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
pQTL_besd <- "/work/database/MR/QTL/pQTL/WholeBlood_UKBPPP_2024/SMR/All"
# work_path <- "/work/run/songzx/MR/Project/MitoDys_SMR_DSZ/3_IVs_Extract"
# coloc_outPath <- "/work/run/songzx/MR/Project/MitoDys_SMR_DSZ/4_IVs_coloc"
# files_id <- "M13_ANKYLOSPON"
# type <- "pQTL"
# coloc_windows <- 100
work_path <- args[1]
coloc_outPath <-  args[2]
files_id <- args[3]
type <- args[4]
coloc_windows <- args[5]
coloc_windows <- as.numeric(coloc_windows)
files_path <- file.path(work_path,files_id,type,"StepII")
coloc_outPath <- file.path(coloc_outPath,files_id,coloc_windows)


dat_pQTL_anno <- read.table("/work/database/MR/QTL/pQTL/WholeBlood_Decode_2021/SMR/5e-8_1000_01/ID_Map_mito.txt",sep = "\t",header = 1,quote="")
dat_eQTL_anno <- read.table("/work/database/MR/QTL/eQTL/eQTLGen_PhaseI/SMR/5e-8_1000_01/ID_Map_mito.txt",sep = "\t",header = 1,quote="")
dat_mQTL_anno <- read.table("/work/database/MR/QTL/mQTL/WholeBlood_2024/EUR/SMR/5e-8_1000_01/ID_Map_mito.txt",sep = "\t",header = 1,quote="")

QTL_anno <- list(
  "mQTL"=dat_mQTL_anno,
  "pQTL"=dat_pQTL_anno,
  "eQTL"=dat_eQTL_anno
)
if(!dir.exists(coloc_outPath)){
  dir.create(coloc_outPath,recursive = T)
}

files <- list.files(files_path,pattern = "*.txt",full.names = T)
QTL_dat <- Reduce(rbind,lapply(files,function(file){
  read.table(file,sep = "\t",header = 1)
}))


if(type == "mQTL"){
  dat_anno <- QTL_anno[[type]] %>% dplyr::select(Gene,GeneBp,Chr)
  Gene2Chr <- lapply(split(dat_anno$Chr,dat_anno$Gene),function(x) unique(x))
  Gene2Start <- lapply(split(dat_anno$GeneBp,dat_anno$Gene),function(x) x-1000*1024)
  Gene2End <- lapply(split(dat_anno$GeneBp,dat_anno$Gene),function(x) x+1000*1024)
  query_res <- QTL_dat %>% dplyr::select(all_of(selCol)) %>%
    dplyr::mutate(start =POS38-coloc_windows*1024,
                  end = POS38+coloc_windows*1024)
  query_res$snp_status <- "trans"
  query_res[apply(query_res,1,function(row) {
    (row["chr"] == Gene2Chr[[row["exposure"]]]) &&(row["POS38"] > Gene2Start[[row["exposure"]]] ) &&
      (row["POS38"] < Gene2End[[row["exposure"]]])
  }),"snp_status"] <- "cis"
  write.table(query_res,sprintf("%s/%s_%d_IVs_coloc.txt",coloc_outPath,type,coloc_windows),row.names = F,quote = F,sep = "\t")
  QTL_leader_SNP_loc <- query_res %>% dplyr::select(id,chr,exposure,SNP) %>%
    dplyr::distinct()
  
  
  chr2ProbeID <- split(QTL_leader_SNP_loc$id,QTL_leader_SNP_loc$chr)
  lapply(names(chr2ProbeID),function(x){
    write.table(unique(chr2ProbeID[[x]]),sprintf("%s/%s_%d_probe.list",coloc_outPath,type,as.numeric(x)),row.names = F,col.names = F,quote = F)
  })
  
  write.table(paste(lapply(names(chr2ProbeID),function(x){
    beqtl <- sprintf("%s/%s_chr%d",mQTL_besd,type,as.numeric(x))
    probe_list <- sprintf("%s/%s_%d_probe.list",coloc_outPath,type,as.numeric(x))
    sprintf("smr --beqtl-summary %s --extract-probe %s --query 1 --out %s_Query_chr%d",beqtl,probe_list,type,as.numeric(x))
  }),sep  = ""),sprintf("%s/Extract_probe_command_%s",coloc_outPath,type),sep = "\t",row.names = F,col.names = F,quote = F)
  
}else if(type == "eQTL"){
  dat_anno <- QTL_anno[[type]] %>% dplyr::select(Gene,ProbeBp,Chr)
  Gene2Chr <- lapply(split(dat_anno$Chr,dat_anno$Gene),function(x) unique(x))
  Gene2Start <- lapply(split(dat_anno$ProbeBp,dat_anno$Gene),function(x) x-1000*1000)
  Gene2End <- lapply(split(dat_anno$ProbeBp,dat_anno$Gene),function(x) x+1000*1000)
  
  query_res <- QTL_dat %>% dplyr::select(all_of(selCol)) %>%
    dplyr::mutate(start =POS38-coloc_windows*1024,
                  end = POS38+coloc_windows*1024)
  query_res$snp_status <- "trans"
  query_res[apply(query_res,1,function(row) {
    (row["chr"] == Gene2Chr[[row["exposure"]]]) &&(row["POS38"] > Gene2Start[[row["exposure"]]] ) &&
      (row["POS38"] < Gene2End[[row["exposure"]]])
  }),"snp_status"] <- "cis"
  write.table(query_res,sprintf("%s/%s_%d_IVs_coloc.txt",coloc_outPath,type,coloc_windows),row.names = F,quote = F,sep = "\t")
  QTL_leader_SNP_loc <- query_res %>% dplyr::select(id,chr,exposure,SNP) %>%
    dplyr::distinct()
  
  chr2ProbeID <- split(QTL_leader_SNP_loc$id,QTL_leader_SNP_loc$chr)
  lapply(names(chr2ProbeID),function(x){
    write.table(unique(chr2ProbeID[[x]]),sprintf("%s/%s_%d_probe.list",coloc_outPath,type,as.numeric(x)),row.names = F,col.names = F,quote = F)
  })
  
  write.table(paste(lapply(names(chr2ProbeID),function(x){
    beqtl <- sprintf(eQTL_besd)
    probe_list <- sprintf("%s/%s_%d_probe.list",coloc_outPath,type,as.numeric(x))
    sprintf("smr --beqtl-summary %s --extract-probe %s --query 1 --out %s_Query_chr%d",beqtl,probe_list,type,as.numeric(x))
  }),sep  = ""),sprintf("%s/Extract_probe_command_%s",coloc_outPath,type),sep = "\t",row.names = F,col.names = F,quote = F)
  
}else if(type == "pQTL"){
  dat_anno <- QTL_anno[[type]] %>% dplyr::select(Gene,ProbeBp,Chr)
  Gene2Chr <- lapply(split(dat_anno$Chr,dat_anno$Gene),function(x) unique(x))
  Gene2Start <- lapply(split(dat_anno$ProbeBp,dat_anno$Gene),function(x) x-1000*1000)
  Gene2End <- lapply(split(dat_anno$ProbeBp,dat_anno$Gene),function(x) x+1000*1000)
  
  query_res <- QTL_dat %>% dplyr::select(all_of(selCol)) %>%
    dplyr::mutate(start =POS38-coloc_windows*1024,
                  end = POS38+coloc_windows*1024)
  query_res$snp_status <- "trans"
  query_res[apply(query_res,1,function(row) {
    (row["chr"] == Gene2Chr[[row["exposure"]]]) &&(row["POS38"] > Gene2Start[[row["exposure"]]] ) &&
      (row["POS38"] < Gene2End[[row["exposure"]]])
  }),"snp_status"] <- "cis"
  write.table(query_res,sprintf("%s/%s_%d_IVs_coloc.txt",coloc_outPath,type,coloc_windows),row.names = F,quote = F,sep = "\t")
  QTL_leader_SNP_loc <- query_res %>% dplyr::select(id,chr,exposure,SNP) %>%
    dplyr::distinct()
  
  probe2SNP <- split(QTL_leader_SNP_loc$SNP,QTL_leader_SNP_loc$id)
  write.table(paste(lapply(names(probe2SNP),function(x){
    beqtl <- sprintf("%s/%s",pQTL_besd,x)
    paste(
      lapply(probe2SNP[[x]],function(snp) sprintf("smr --beqtl-summary %s --snp %s --query 1  --snp-wind %d --out %s_Query_%s_%s",beqtl,snp,coloc_windows,type,snp,x))
    ,collapse = "\n")
  }),sep  = ""),sprintf("%s/Extract_probe_command_%s",coloc_outPath,type),sep = "\t",row.names = F,col.names = F,quote = F)
  
}

