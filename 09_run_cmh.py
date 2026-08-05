# just PTV as example

hailctl dataproc submit [machine] SCHEMA_counts_export.py \
--mt [MT PATH] \
--vep_ht [VEP PATH] \
--cons_cat ptv \
--ac 15 \
--out cmh/ \
--file_prefix [FILE NAME] \
--ac_file [ALLELE COUNT PATH]

hailctl dataproc submit js run_cmh_updated_with_cis.py \
--mt [MT FROM COUNTS EXPORT] \
--manifest [META PATH] \
--counts [COUNTS DATA PATH] \
--ptv_only \
--out [RESULTS PATH]
