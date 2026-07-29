import hail as hl
hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")
hl.default_reference('GRCh38')
from gnomad.utils.vep import vep_or_lookup_vep

MT = 'schema2_scz_case_control_11-04-2025.mt'

ht = hl.read_matrix_table(MT).rows()

ht = ht.repartition(75000)
ht.write('schema2_scz_case_control_vars_to_vep_11-05-2025.ht', overwrite=True)

print('read back in')
ht = hl.read_table('schema2_scz_case_control_vars_to_vep_11-05-2025.ht')

print(f"Generating VEP annotations...")
ht_vep = vep_or_lookup_vep(ht, vep_version=95)


# Write VEP annotations to HT output
OUT = 'schema2_scz_case_control_vep_11-05-2025.ht'

print(f"Writing VEP annotations to {OUT}...")
ht_vep.write(OUT, overwrite = True)

