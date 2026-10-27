Revision calculations for the inflammatory arthritis manuscript

Run with Python 3: python revision_calculations.py
No third-party packages or network access are required. Results are written to calculated/.

Scope
1. Reweight the original precise FinnGen PPH4 values for p12 = 1e-6, 5e-6, 1e-5 and 5e-5, holding p1 = p2 = 1e-4 and the regional evidence fixed.
2. Calculate 95% CIs, QTL F statistics and effect-ratio discrepancies for the two supplied 47-row VARS2 outputs.

This code does not rerun the original GWAS, SMR/HEIDI, regional colocalization or disease-wide FDR. It does not perform conditional analysis or recover missing annotations. The 47-row files are retained outputs; no complete regional or disease-wide testing family is assumed.

Input provenance
finnGen_coloc_input.csv preserves precise PPH4 and original q/HEIDI information underlying Supplementary Tables 2, 4 and 6. baseline_primary means reported q < 0.05, P_HEIDI > 0.05 and PPH4 > 0.70.
The two .smr files are unchanged copies supplied with the 27 October 2026 analysis.

Posterior update
For h = PPH4 at baseline p12 = 1e-5 and r = new_p12 / 1e-5, updated PPH4 = r*h / (1-h+r*h). Only the shared-variant prior changes. Source: https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html

Outcome verification
Neale Lab 2018 both-sexes field 20002 code 1477 is self-reported psoriatic arthropathy, with 712 cases and 360429 controls. The rs6920760 record has beta 0.000997816, SE 0.000129117 and P 1.09514e-14, exactly matching the supplied SMR output. In the public GRCh37 file it is 6:30874797:C:T; the SMR file uses different genomic coordinates. Match uses allele-oriented statistics, not an assumption that the builds agree.
Metadata: https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.both_sexes.tsv.bgz
GWAS: https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1477.gwas.imputed_v3.both_sexes.tsv.bgz
