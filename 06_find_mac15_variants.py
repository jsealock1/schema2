import hail as hl

mt = hl.read_matrix_table(MT)

vqc = hl.variant_qc(mt).rows()
print('filter to mac')
mac = vqc.annotate(MAC = vqc.variant_qc.AC[1])
mac = mac.select('MAC')

schema1_totals = hl.read_table(SCHEMA1)

## combine ipysch and callset macs
callset_vars = mac.select()
schema1_vars = schema1_totals.select()

schema1_vars_unique = schema1_vars.filter(hl.is_defined(callset_vars[schema1_vars.key]), keep=False)

ht = callset_vars.union(schema1_vars_unique)
ht = ht.distinct()

ht = ht.annotate(callset_mac = mac[ht.key].MAC, ipsych_uk10k_mac = schema1_totals[ht.key].MAC)
ht = ht.annotate(callset_mac = hl.coalesce(ht.callset_mac, 0), ipsych_uk10k_mac = hl.coalesce(ht.ipsych_uk10k_mac, 0))
ht = ht.annotate(MAC = ht.callset_mac + ht.ipsych_uk10k_mac)

print('repart and write all')
ht = ht.filter(ht.MAC > 0)

print('write')
ht.write(OUT_ALL, overwrite=True)
