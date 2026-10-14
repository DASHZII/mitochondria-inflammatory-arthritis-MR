import polars as pl
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

def Extract_coloc(PosDict,Query_res,probe_query,bp_query,start="start",end ="end"):
    dat_res = {}
    empty_index = []
    for k,v in  PosDict.items():
        for patch in v:
            tmp_res = Query_res.filter([
                pl.col(probe_query) == k,
                pl.col(bp_query) > patch[start], 
                pl.col(bp_query) < patch[end]
                ]  
            ).with_columns(
                pl.lit(patch["SNP"]).alias("Source")
            )
            if tmp_res.shape[0] == 0:
                empty_index.append(k)
            else:
                dat_res["_".join([k,patch["SNP"]])] = tmp_res
    return(dat_res,empty_index)

def Coordinary_Trans(chr,pos,lo):
    # 使用convert_coordinate进行坐标转换  
    converted = lo.convert_coordinate(f"chr{chr}", pos, pos)
    if converted:  
        loc = converted[0][1] 
    else:  
       loc = None
    return(loc)

def Merge_QTL_GWAS(datQTL,datGWAS,samplesize,probe,SNP):

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
        pl.lit(probe).alias("probe"),
        pl.lit(SNP).alias("source")
    )
    return(merge_res)
    
samplesizeD = {
    "eQTL":31864,
    "mQTL":3701
}

versionD={
    "eQTL":"GRCh37",
    "mQTL":"GRCh37",
    "pQTL":"GRCh38"
}
pQTL_map = pl.read_csv("/work/database/MR/QTL/pQTL/WholeBlood_Decode_2021/Decode_FileId_map_full_sample.txt",separator="\t")
id2N = dict(zip(pQTL_map["SeqId"],pQTL_map["n"]))
chainfile = '/work/database/MR/Reference/Chain/hg19ToHg38.over.chain'


gwas_path = sys.argv[1]
coloc_Path = sys.argv[2]
file_id = sys.argv[3]
type = sys.argv[4]
coloc_windows=sys.argv[5]

version = versionD[type]



inputPath = Path(coloc_Path).joinpath(file_id).joinpath(coloc_windows)

dat_sig = pl.read_csv(f"{inputPath}/{type}_{coloc_windows}_IVs_coloc.txt",separator="\t")
Query_res = pl.concat([pl.read_csv(file,separator="\t",schema_overrides=col_default_type,ignore_errors=True) for file in Path(inputPath).glob("_".join([type,"Query*"]))])
disease_gwas = pl.read_csv(f"{gwas_path}/{file_id}_gwas.ma",separator="\t")


Query_res = Query_res.filter(
    ~Query_res.is_duplicated()
)

if version == "GRCh37":
    # 初始化liftOver对象  
    lo = pyliftover.LiftOver(chainfile)  
    # 将37的坐标系转化为38的
    Query_res = Query_res.with_columns(
        Query_res.map_rows(lambda row:Coordinary_Trans(row[1],row[2],lo))
    ).drop(["BP"]).rename({"map":"BP"}).drop_nulls()

mapping = {}
for row in dat_sig.iter_rows():
    id_value = row[4]
    if id_value not in mapping:
        mapping[id_value] = []
    mapping[id_value].append({
        "SNP": row[0],
        "start": row[6],
        "end": row[7]
    })

dat_res,empty_index = Extract_coloc(mapping,Query_res,"Probe","BP")

disease_gwas = disease_gwas.with_columns(
    pl.when(pl.col("freq") > 0.5).\
    then(1-pl.col("freq")).\
    otherwise(pl.col("freq")).alias("MAF")
).select(["SNP","MAF","b","se","p","n"])

if type == "pQTL":
    dat_coloc = pl.concat([Merge_QTL_GWAS(v,disease_gwas,id2N["_".join((k.split("_")[0:2]))],"_".join((k.split("_")[0:2])),k.split("_")[2]) for k,v in dat_res.items()])
else:
    dat_coloc = pl.concat([Merge_QTL_GWAS(v,disease_gwas,samplesizeD[type],k.split("_")[0],k.split("_")[1]) for k,v in dat_res.items()])
dat_coloc.write_csv(f"{coloc_Path}/{file_id}/coloc_{type}_{coloc_windows}.txt",separator="\t")