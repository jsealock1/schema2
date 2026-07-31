import hail as hl

## DD/ID
genes = hl.import_table(GENES, key='gene')
dd_genes = genes.filter(genes.phenotype=='DD')

mt_plof = hl.read_matrix_table(MT_PLOF)
mt_plof_dd = mt_plof.filter_rows(hl.is_defined(dd_genes[mt_plof.gene_symbol]))
out_plof_dd = (mt_plof_dd.annotate_cols(total_plof_os = hl.agg.sum(mt_plof_dd.agg)).cols())

mt_mis = hl.read_matrix_table(MT_MIS)
mt_mis_dd = mt_mis.filter_rows(hl.is_defined(dd_genes[mt_mis.gene_symbol]))
out_mis_dd = (mt_mis_dd.annotate_cols(total_misrank93 = hl.agg.sum(mt_mis_dd.agg)).cols())

mt_syn = hl.read_matrix_table(MT_SYN)
mt_syn_dd = mt_syn.filter_rows(hl.is_defined(dd_genes[mt_syn.gene_symbol]))
out_syn_dd = (mt_syn_dd.annotate_cols(total_syn = hl.agg.sum(mt_syn_dd.agg)).cols())

out_dd = out_plof_dd.annotate(total_misrank93 = out_mis_dd[out_plof_dd.s].total_misrank93, 
                                total_syn = out_syn_dd[out_plof_dd.s].total_syn)

print('out dd')
out_dd.export(output=OUT_DD)



## ASD
mt_plof = hl.read_matrix_table(MT_PLOF)
mt_plof_asd = mt_plof.filter_rows(hl.is_defined(asd_genes[mt_plof.gene_symbol]))
out_plof_asd = (mt_plof_asd.annotate_cols(total_plof_os = hl.agg.sum(mt_plof_asd.agg)).cols())

mt_mis = hl.read_matrix_table(MT_MIS)
mt_mis_asd = mt_mis.filter_rows(hl.is_defined(asd_genes[mt_mis.gene_symbol]))
out_mis_asd = (mt_mis_asd.annotate_cols(total_misrank93 = hl.agg.sum(mt_mis_asd.agg)).cols())

mt_syn = hl.read_matrix_table(MT_SYN)
mt_syn_asd = mt_syn.filter_rows(hl.is_defined(asd_genes[mt_syn.gene_symbol]))
out_syn_asd = (mt_syn_asd.annotate_cols(total_syn = hl.agg.sum(mt_syn_asd.agg)).cols())

out_asd = out_plof_asd.annotate(total_misrank93 = out_mis_asd[out_plof_asd.s].total_misrank93, 
                                total_syn = out_syn_asd[out_plof_asd.s].total_syn)

print('out asd')
out_asd.export(output=OUT_ASD)


