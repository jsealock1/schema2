import hail as hl

hl.default_reference('GRCh38')
from gnomad.utils.vep import vep_or_lookup_vep

ht = hl.read_matrix_table(MT).rows()

ht = ht.repartition(75000)

print(f"Generating VEP annotations...")
ht_vep = vep_or_lookup_vep(ht, vep_version=95)


# Write VEP annotations to HT output

print(f"Writing VEP annotations to {OUT}...")
ht_vep.write(OUT, overwrite = True)

