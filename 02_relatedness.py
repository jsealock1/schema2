import hail as hl

## relatedness
mt = hl.read_matrix_table(MT)

mt = hl.variant_qc(mt)
mt = mt.filter_rows(hl.len(mt.alleles) == 2)
mt = mt.filter_rows(((mt.variant_qc.AF[0] > 0.01) & (mt.variant_qc.AF[1] > 0.01)) & ((mt.variant_qc.AF[0] < 0.99) & (mt.variant_qc.AF[1] < 0.99)))
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.call_rate>0.9)

pca_snps = hl.read_table(PCA_SNPS)
mt2 = mt.filter_rows(hl.is_defined(pca_snps[mt.row_key]), keep=True)

dataset = mt2
rel = hl.pc_relate(dataset.GT, 0.01, k=10, statistics='kin', min_kinship=0.05) 

print('save relatedness')
rel.write(REL_OUT, overwrite=True)

samples = dataset.cols()

pairs = rel.filter(rel['kin'] > 0.125)
related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
print(related_samples_to_remove.count())






