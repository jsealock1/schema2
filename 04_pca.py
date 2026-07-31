
import hail as hl


from gnomad.sample_qc import *
from gnomad.sample_qc.ancestry import *



##### calclated PCS
MT = 'schema2_scz_case_control_gnomadv4_snps_11-18-2025_call_rate_repart.mt'

import hail as hl

mt = hl.read_matrix_table(MT)

print('calc pcs')
_, scores, _ = hl.hwe_normalized_pca(mt.GT)

OUT = 'schema2_scz_case_control_gnomadv4_snps_calculated_pcs_all.ht'
OUT_TSV = 'schema2_scz_case_control_gnomadv4_snps_calculated_pcs_all.tsv.bgz'

print('write')
scores.write(OUT)
scores.flatten().export(OUT_TSV)

