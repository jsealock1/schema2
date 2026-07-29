# just PTV as example
hailctl dataproc start js --num-workers 2 --max-idle=45m --master-machine-type n1-highmem-16 --worker-machine-type n1-highmem-16 --autoscaling-policy=max-50 --public-ip-address

hailctl dataproc submit js SCHEMA_counts_export.py \
--mt schema2_scz_case_control_unrelated_lcr_geno_dp_filtered_11-20-2025.mt \
--vep_ht schema2_scz_case_control_vep_annotated_with_cmh_annos_11-10-2025.ht \
--cons_cat ptv \
--ac 15 \
--out cmh/ \
--file_prefix 2025_11_21_ALL_SCHEMA_PTV_MAC15_GQ25_DPMAX_NO_UNK \
--ac_file schema2_scz_cases_controls_all_vars_mac_11-20-2025.ht

hailctl dataproc submit js run_cmh_updated_with_cis.py \
--mt 2025_11_21_ALL_SCHEMA_PTV_MAC15_GQ25_DPMAX_NO_UNK_ptv_AC15_counts.mt \
--manifest schema2_pc_matched_01-07-2026.tsv \
--counts ipsych_uk10k_trios_all_variants_annotated_MAC15_raw_misrank_11-10-2025.ht \
--ptv_only \
--out 2026_06_12_ALL_SCHEMA_PTV_MAC15_GQ25_ptv_AC15_CMH_pc_matched_fishers.tsv
