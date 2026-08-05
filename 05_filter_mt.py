import hail as hl

## filter to unrelated

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


mt.write(OUT2, overwrite=True)
