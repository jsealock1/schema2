import hail as hl
hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")

## filter to unrelated
MT = 'gs://schema_jsealock/bge_wave3/schema2_add_dragen/schema2_scz_case_control_11-11-2025.mt'
RELATED_SAMPLES = 'gs://schema_jsealock/schema2/november_2025/schema2/data_king_related_samples_to_remove.tsv'
KEY = 'gs://schema_jsealock/schema2/november_2025/schema2/king_relatedness_sample_key.tsv'

mt = hl.read_matrix_table(MT)
related_samples_to_remove = hl.import_table(RELATED_SAMPLES, key='node')
key = hl.import_table(KEY, key='s2')
related_samples_to_remove = related_samples_to_remove.annotate(s = key[related_samples_to_remove.node].s)


related_samples_to_remove = related_samples_to_remove.key_by('s')
mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.s]), keep=False)

## 

## remove indels, lcrs, bad genotype calls

LCR = 'gs://gcp-public-data--gnomad/resources/grch38/lcr_intervals/LCRFromHengHg38.ht'

print('filter lcrs')
lcr = hl.read_table(LCR)

intervals_to_exclude = lcr.interval.collect()
mt = hl.filter_intervals(mt, intervals_to_exclude, keep=False)

print(' remove indels > 10bp')

mt = mt.filter_rows(
    hl.case()
    .when(mt.alleles[0].length() != mt.alleles[1].length(),  # Check for indels
          hl.abs(mt.alleles[0].length() - mt.alleles[1].length()) <= 10)
    .default(True)  # Keep other variants (e.g., SNPs)
)



print('filter to AB >= 0.25 and GQ >= 25 for rvas')
mt = mt.filter_entries(
        hl.is_defined(mt.GT) &
        (
            (mt.GT.is_het() & (((mt.AD[1] / mt.DP) < 0.25) | (mt.GQ < 25))) |
            (mt.GT.is_hom_var() & (((mt.AD[1] / mt.DP) < 0.8) | (mt.GQ < 25)))
    ),
    keep=False
)


mt = mt.filter_entries(((mt.capture == 'Nextera') & (mt.DP < 12)), keep=False)


mt = mt.filter_entries(
    hl.case()
      .when(
          (mt.callset == 'dragen') |
          (mt.callset == 'bge') |
          (mt.callset == 'gnomad_genomes'),
          mt.DP > 80
      )
      .when(
          (mt.callset == 'cardiff') |
          (mt.callset == 'gnomad_exomes'),
          mt.DP > 200
      )
      .default(False),
    keep=False
)

OUT2 = 'schema2_scz_case_control_unrelated_lcr_geno_dp_filtered_11-19-2025.mt'
mt.write(OUT2, overwrite=True)


## pc matched, no unknown cohort
MT = 'schema2_scz_case_control_unrelated_lcr_geno_dp_filtered_11-19-2025.mt'
META = 'schema2_pc_matched.tsv'

import hail as hl
hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")

mt = hl.read_matrix_table(MT)
meta = hl.import_table(META, key='s')

mt = mt.filter_cols(hl.is_defined(meta[mt.col_key]))
mt.count_cols()

OUT3 = 'schema2_scz_case_control_unrelated_lcr_geno_dp_filtered_11-20-2025.mt'
mt.write(OUT3, overwrite=True)
