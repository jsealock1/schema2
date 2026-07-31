## annotate vep with mis mean rank

import hail as hl

vep = hl.read_table(VEP)
vep_dn = hl.read_table(DN_VEP)

vep_dn2 = vep_dn.filter(hl.is_defined(vep[vep_dn.key]), keep=False)

ht = hl.read_table(HT)
ht = ht.order_by(ht.mean_rank)
ht = ht.add_index(name="rank")  
total_rows = ht.count()
ht = ht.annotate(mean_rank_percentile=(ht.rank / total_rows) * 100)
ht = ht.key_by('locus','alleles')

vep = hl.read_table(VEP)


vep = vep.annotate(mis_mean_rank = ht[vep.key].mean_rank_percentile)

vep = vep.drop('consequence_category','vep')
vep = vep.annotate(consequence_category = hl.if_else(vep.isPTV, 'pLoF', 
                                            hl.if_else(vep.isLC, 'LC',
                                            hl.if_else(vep.isSYN, 'synonymous',
                                            hl.if_else(vep.isMIS, 'missense', 'other')))))

vep.write(OUT, overwrite=True)
