import hail as hl

print('load data')
dragen = hl.read_matrix_table(DRAGEN)

## filter on call rate
dragen = dragen.select_rows()
dragen = hl.variant_qc(dragen)
dragen = dragen.filter_rows(dragen.variant_qc.call_rate >= 0.95)


cardiff_mt = hl.read_matrix_table(CARDIFF)
gnomad_mt = hl.read_matrix_table(GNOMAD)

bge_mt = hl.read_matrix_table(BGE)
bge_mt = bge_mt.select_rows()
bge_mt = hl.variant_qc(bge_mt)
bge_mt = bge_mt.filter_rows(bge_mt.variant_qc.call_rate >= 0.95)

wgspd = hl.read_matrix_table(WGSPD)

print("union")
bge_mt = bge_mt.select_entries('GQ','DP','GT','AD').select_cols().select_rows()
cardiff_mt = cardiff_mt.select_entries('GQ','DP','GT','AD').select_cols().select_rows()
gnomad_mt = gnomad_mt.select_rows().select_entries('GQ','DP','GT','AD').select_cols().select_rows()
wgspd = wgspd.select_entries('GQ','DP','GT','AD').select_cols().select_rows()
dragen = dragen.select_entries('GQ','DP','GT','AD').select_cols().select_rows()

mt1 = bge_mt.union_cols(cardiff_mt, row_join_type='outer') # bge cardiff
mt2 = mt1.union_cols(gnomad_mt, row_join_type='outer') # bge cardiff gnomad ex
mt3 = mt2.union_cols(wgspd, row_join_type='outer') # bge, cardiff, gnomad ex, gnomad wgs
mt = mt3.union_cols(dragen, row_join_type='outer') # bge, cardiff, gnomad ex, gnomad wgs, dragen


## remove related samples

related_samples_to_remove = hl.import_table(RELATED_SAMPLES, key='node')
key = hl.import_table(KEY, key='s2')
related_samples_to_remove = related_samples_to_remove.annotate(s = key[related_samples_to_remove.node].s)

related_mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.s]))
related_mt = related_mt.cols()

mt = mt.filter_cols(~hl.is_defined(related_mt[mt.col_key]))

## filter to high qual scz case control
meta = hl.import_table(META, key='s')

meta = meta.filter(meta.high_quality == 'true')
meta = meta.filter((meta.status == 'control') | (meta.status=='schizophrenia'))

mt = mt.filter_cols(hl.is_defined(meta[mt.s]))
mt = mt.annotate_cols(**meta[mt.s])

print("write")
mt.write(OUT, overwrite=True)
