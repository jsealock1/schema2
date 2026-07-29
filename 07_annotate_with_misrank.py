## annotate vep with mis mean rank

import hail as hl
hl.init(gcs_requester_pays_configuration = 'daly-neale-sczmeta', default_reference = 'GRCh38', tmp_dir='gs://schema_jsealock/tmp/')
# make annotations usable for cmh
    # ht = ht.annotate(consequence = hl.if_else((ht.consequence_category == "pLoF") | (ht.other_splice), 'ptv',
    #                             hl.if_else(ht.consequence_category == "LC", 'ptv_lc',
    #                             hl.if_else((ht.consequence_category == "synonymous") & (ht.other_splice==False), 'synonymous',
    #                             hl.if_else(((ht.mis_mean_rank >= 93) & (ht.other_splice==False)),  'mis_mean_rank93', 
    #                             'none')))))

VEP = 'gs://schema_jsealock/schema2/november_2025/vep/schema2_scz_case_control_vep_annotated_11-05-2025.ht'
DN_VEP = 'gs://schema_jsealock/de_novo_analysis/schema1_de_novo_variants_grch38_vep_annotated.ht'

vep = hl.read_table(VEP)
vep_dn = hl.read_table(DN_VEP)

vep_dn2 = vep_dn.filter(hl.is_defined(vep[vep_dn.key]), keep=False)


HT = 'gs://schema_jsealock/schema2/add_dragen/jan_2025/cmh/missense_ranking/2025-03-24_missense-annotations_mean-rank_all.ht'
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

OUT = 'gs://schema_jsealock/schema2/november_2025/vep/schema2_scz_case_control_vep_annotated_with_cmh_annos_11-10-2025.ht'
vep.write(OUT, overwrite=True)