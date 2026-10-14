import yaml
import json
from pathlib import Path
import polars as pl
import pandas as pd
import numpy as np
from functools import reduce
import sys

from MR_Py import *

configF = None
sup_diction = Path("/work/database/MR/Reference/sup")
if configF is None:
    configF = sup_diction/"config.yaml"
with open(configF,"r",encoding="UTF-8") as f:
    config = yaml.load(f.read(),Loader=yaml.FullLoader)

with open(sup_diction/'output.json', 'r',encoding="utf-8") as file:  
    Protein_default_col = json.load(file)


versionD={
    "eQTL":"GRCh37",
    "mQTL":"GRCh37",
    "pQTL":"GRCh38"
}
samplesizeD = {
    "eQTL":31864,
    "mQTL":3701
}
pQTL_map = pl.read_csv("/work/database/MR/QTL/pQTL/WholeBlood_UKBPPP_2024/SMR/ID_Map_mito2.txt",separator="\t")
id2N = dict(zip(pQTL_map["Probe"],pQTL_map["n"]))

Inpath = sys.argv[1]
qtl_type = sys.argv[2]


if ~os.path.exists(f"{Inpath}/{qtl_type}/StepII"):
    os.makedirs(f"{Inpath}/{qtl_type}/StepII",exist_ok=True)
p_threshold = 5e-8
F_threshold = 10

inputpath = Path(f"{Inpath}/{qtl_type}/StepI")

files = [file for file in inputpath.glob(pattern="*.txt")]

se_col = "SE" if  qtl_type != "pQTL" else "se"
empty_list = []
for file in files:
    filename =  os.path.basename(file).replace(".txt","")
    phenotype_id =filename.split("_")[0]
    file_id = "_".join(filename.split("_")[1:3])

    samplesize_num = samplesizeD[qtl_type] if qtl_type != "pQTL" else id2N[filename]

    exposure_data = ExposureProcess(infile=file,build = versionD[qtl_type],p_threshold = p_threshold, f_threshold=F_threshold, config=config,
                            snp_col="SNP", chr_col="Chr", pos_col="BP", other_allele_col="A2",
                            effect_allele_col="A1", eaf_col="Freq", maf_col=None, beta_col="b",
                            se_col=se_col, pval_col="p", sample_size=None, or_col=None,
                            sample_size_num=samplesize_num,pattern = None,id_=file_id,phenotype=phenotype_id)
    dat_exp,dat_res = exposure_data.exposure()
    if dat_exp is not None:
        dat_exp.write_csv(Path(f"{Inpath}/{qtl_type}/StepII/{filename}_IVs.txt"),separator="\t")
    else:
        empty_list.append(filename)

