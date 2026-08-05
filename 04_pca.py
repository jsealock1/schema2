
import hail as hl


from gnomad.sample_qc import *
from gnomad.sample_qc.ancestry import *

mt = hl.read_matrix_table(MT)

print('calc pcs')
_, scores, _ = hl.hwe_normalized_pca(mt.GT)

print('write')
scores.write(OUT)
scores.flatten().export(OUT_TSV)

