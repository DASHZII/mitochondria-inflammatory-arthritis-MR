import polars as pl
import glob
from pathlib import Path
import pyliftover  
import sys

col_default_type = {
    "SNP":pl.String,
    "BP":pl.Int64,
    "Freq":pl.Float64,
    "b":pl.Float64,
    "SE":pl.Float64,
    "p":pl.Float64
}

def Extract_coloc(dat_meta,Query_res):
    start_pos = dict(zip(dat_meta["probeID"],dat_meta["start"]))
    end_pos = dict(zip(dat_meta["probeID"],dat_meta["end"]))
    dat_res = {}
    empty_index = []
    for probe in dat_meta["probeID"].to_list():
        tmp_res = Query_res.filter([
            pl.col("Probe") == probe,
            pl.col("BP") > start_pos[probe], 
            pl.col("BP") < end_pos[probe]
            ]  
        )
        if tmp_res.shape[0] == 0:
            empty_index.append(probe)
        else:
            dat_res[probe] = tmp_res
    return(dat_res,empty_index)

def Coordinary_Trans(chr,pos,lo):
    # 使用convert_coordinate进行坐标转换  
    converted = lo.convert_coordinate(f"chr{chr}", pos, pos)
    if converted:  
        loc = converted[0][1] 
    else:  
       loc = None
    return(loc)

def Merge_QTL_GWAS(datQTL,datGWAS,samplesize,probe):
    tmp = datQTL.with_columns(
        pl.lit(samplesize).alias("n"),
        pl.when(pl.col("Freq") > 0.5).\
        then(1-pl.col("Freq")).\
        otherwise(pl.col("Freq")).alias("MAF")
    ).select(["SNP","BP","p","MAF","b","SE","n"]).\
    rename({"BP":"bp","SE":"se"})
    merge_res = tmp.rename({x:"_".join([x,"QTL"]) for x in tmp.columns if x not in ["SNP","bp"]}).\
        join(datGWAS.rename({x:"_".join([x,"GWAS"]) for x in datGWAS.columns if x != "SNP"}),on= "SNP")
    merge_res = merge_res.with_columns(
        pl.lit(probe).alias("probe")
    )
    return(merge_res)
    
samplesizeD = {
    "eQTL":31864,
    "mQTL":3701
}
pQTL_map = pl.read_csv("/work/database/MR/QTL/pQTL/WholeBlood_Decode_2021/Decode_FileId_map_full_sample.txt",separator="\t")
id2N = dict(zip(pQTL_map["SeqId"],pQTL_map["n"]))


# version = "GRCh37"
# type = "pQTL"
# file_id = "M45"
# gwas_path = "/work/run/songzx/MR/Project/MitoDys_SMR_DSZ/raw_data/gwas-summary"
# coloc_Path = "/work/run/songzx/MR/Project/MitoDys_SMR_DSZ/2_coloc_res"
# chainfile = '/work/database/MR/Reference/Chain/hg19ToHg38.over.chain'
gwas_path = sys.argv[1]
coloc_Path = sys.argv[2]
file_id = sys.argv[3]
type = sys.argv[4]
version = sys.argv[5]
chainfile = '/work/database/MR/Reference/Chain/hg19ToHg38.over.chain'

# 初始化liftOver对象  
lo = pyliftover.LiftOver(chainfile)  
inputPath = Path(coloc_Path).joinpath(file_id)

dat_sig = pl.read_csv(f"{coloc_Path}/{file_id}/{type}_sig_meta.txt",separator="\t")
Query_res = pl.concat([pl.read_csv(file,separator="\t",schema_overrides=col_default_type,ignore_errors=True) for file in Path(inputPath).glob("_".join([type,"Query*"]))])
Query_res = Query_res.drop_nulls()
if version == "GRCh37":
    # 将37的坐标系转化为38的
    Query_res = Query_res.with_columns(
        Query_res.map_rows(lambda row:Coordinary_Trans(row[1],row[2],lo))
    ).drop(["BP"]).rename({"map":"BP"})

dat_res,empty_index = Extract_coloc(dat_sig,Query_res)
disease_gwas = pl.read_csv(f"{gwas_path}/{file_id}_gwas.ma",separator="\t")
disease_gwas = disease_gwas.with_columns(
    pl.when(pl.col("freq") > 0.5).\
    then(1-pl.col("freq")).\
    otherwise(pl.col("freq")).alias("MAF")
).select(["SNP","MAF","b","se","p","n"])

if type == "pQTL":
    dat_coloc = pl.concat([Merge_QTL_GWAS(v,disease_gwas,id2N[k],k) for k,v in dat_res.items()])
else:
    dat_coloc = pl.concat([Merge_QTL_GWAS(v,disease_gwas,samplesizeD[type],k) for k,v in dat_res.items()])
dat_coloc.write_csv(f"{coloc_Path}/{file_id}/coloc_{type}.txt",separator="\t")
