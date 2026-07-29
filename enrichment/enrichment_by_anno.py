import hail as hl
hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")
hl._set_flags(use_new_shuffle='1')

mt = hl.read_matrix_table(MT_IN)

# annotate with vep
vep = hl.read_table(HT_ANNOT)
vep = vep.filter(vep.fail_50_bp_rule_and_gerp_dist, keep=False)


mt = mt.annotate_rows(**vep[mt.row_key])

print('use pli instead')
loeuf = hl.read_table(LOEUF)
pli = loeuf.filter(loeuf.lof.pLI > 0.9)
pli = pli.key_by('gene')

mt2 = mt.filter_rows(hl.is_defined(pli[mt.gene_symbol]))

mac = hl.read_table(MAC)
mac = mac.filter(mac.MAC == 1)

mt3 = mt2.filter_rows(hl.is_defined(mac[mt2.row_key]))


print("find total n ac1 per cat")
out1 = mt3.annotate_cols(
    total_plof_os = hl.agg.sum(mt3.isPTV & mt3.GT.is_non_ref()), 
    total_misrank93 = hl.agg.sum((mt3.mis_mean_rank >= 93) & mt3.GT.is_non_ref()), 
    total_syn = hl.agg.sum(mt3.isSYN & mt3.GT.is_non_ref()) 
    ).cols()


print("export")
out1.flatten().export(output=OUT_PLI1)


##  mac 2-5
mac = hl.read_table(MAC)
mac = mac.filter((mac.MAC > 1) & (mac.MAC < 6)) # 2-5

mt2 = hl.read_matrix_table(MT2)
mt3 = mt2.filter_rows(hl.is_defined(mac[mt2.row_key]))

mt3 = mt2.filter_rows(hl.is_defined(mac[mt2.row_key]))


print("find total n ac2-5 per cat")
out2 = mt3.annotate_cols(
    total_plof_os = hl.agg.sum(mt3.isPTV & mt3.GT.is_non_ref()), 
    total_misrank93 = hl.agg.sum((mt3.mis_mean_rank >= 93) & mt3.GT.is_non_ref()), 
    total_syn = hl.agg.sum(mt3.isSYN & mt3.GT.is_non_ref()) 
    ).cols()


print("export")
out2.flatten().export(output=OUT_PLI2)


# ## MAC 6-10
mac = hl.read_table(MAC)
mac = mac.filter((mac.MAC > 5) & (mac.MAC < 11)) # 6-10
mt3 = mt2.filter_rows(hl.is_defined(mac[mt2.row_key]))


print("find total n ac5-10 per cat")
out3 = mt3.annotate_cols(
    total_plof_os = hl.agg.sum(mt3.isPTV & mt3.GT.is_non_ref()), 
    total_misrank93 = hl.agg.sum((mt3.mis_mean_rank >= 93) & mt3.GT.is_non_ref()), 
    total_syn = hl.agg.sum(mt3.isSYN & mt3.GT.is_non_ref()) 
    ).cols()


print("export")
out3.flatten().export(output=OUT_PLI3)


# ## MAC 10-15
mac = hl.read_table(MAC)
mac = mac.filter((mac.MAC > 10) & (mac.MAC < 16)) # 11-15
mt3 = mt2.filter_rows(hl.is_defined(mac[mt2.row_key]))


print("find total n ac11-15 per cat")
out4 = mt3.annotate_cols(
    total_plof_os = hl.agg.sum(mt3.isPTV & mt3.GT.is_non_ref()), 
    total_misrank93 = hl.agg.sum((mt3.mis_mean_rank >= 93) & mt3.GT.is_non_ref()), 
    total_syn = hl.agg.sum(mt3.isSYN & mt3.GT.is_non_ref()) 
    ).cols()



print("export")
out4.flatten().export(output=OUT_PLI4)





# get total syn

import hail as hl
hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")
hl._set_flags(use_new_shuffle='1')

mt = hl.read_matrix_table(MT_IN)
total_singleton = hl.sample_qc(mt).cols()
total_singleton = total_singleton.select(total_singleton.sample_qc.n_singleton)
total_singleton = total_singleton.rename({'n_singleton' : 'n_total_singleton'})

total_singleton.export('counts_of_total_singletons_11_20_2025.tsv.bgz')
