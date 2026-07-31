
import hail as hl

ht_in = hl.read_table(HT)
ht_in = ht_in.filter(ht_in.consequence=='none', keep=False)

mac_in = hl.read_table(MAC)
mac_in = mac_in.filter(mac_in.MAC > 0)

vep = hl.read_table(VEP)



# use mis mean rank
ht_in = ht_in.annotate(mis_mean_rank = vep[ht_in.key].mis_mean_rank, isOS = vep[ht_in.key].isOS)
ht_in = ht_in.annotate(consequence_category = hl.if_else(ht_in.consequence == 'ptv', 'ptv',
                                                hl.if_else(ht_in.isOS, 'ptv', 
                                                    hl.if_else(ht_in.mis_mean_rank >= 93, 'mis_mean_rank93',
                                                        hl.if_else(ht_in.consequence == 'synonymous', 'synonymous', 'none')))))

ht_in = ht_in.filter(ht_in.consequence == 'none', keep=False)


def make_rvas_file(ht, mac_in, mac):
    print(mac)
    mac_vars = mac_in.filter(mac_in.MAC <= mac)
    ht = ht.filter(hl.is_defined(mac_vars[ht.key]))
    # filter by dataset
    taiwan_trios = ht.filter(ht.dataset=='taiwanese_trios')
    bulg_trios = ht.filter(ht.dataset=='bulgarian_trios')
    uk10k = ht.filter(ht.dataset=='uk10k_interval')
    ipsych = ht.filter(ht.dataset=='danish_case_control')
    # in each group, add cohort, platform, pop, status
    taiwan_trios = taiwan_trios.annotate(pop = 'EAS', cohort='taiwan_trios', capture='Agilent')
    bulg_trios = bulg_trios.annotate(pop = 'EUR', cohort='bulgarian_trios', capture='Agilent')
    ipsych = ipsych.filter(ipsych.group=='EUR-N_non-nextera')
    ipsych = ipsych.annotate(pop = 'EUR', cohort='ipsych', capture='Agilent')
    uk10k = uk10k.filter((uk10k.group=='EUR_non-nextera') | (uk10k.group=='FIN_non-nextera')) ## remove non eur samples
    uk10k = uk10k.annotate(cohort='uk10k', capture='Agilent')
    uk10k = uk10k.annotate(pop = hl.if_else(uk10k.group=='EUR_non-nextera', 'EUR','FIN'))
    trios = taiwan_trios.union(bulg_trios)
    trios_ipsych = trios.union(ipsych)
    ht = trios.union(uk10k)
    # combine and add callset=counts, status as string
    ht = ht.annotate(status = hl.if_else(ht.SCZ, 'schizophrenia', 'control'))
    ht = ht.group_by(ht.gene_symbol, ht.consequence_category, ht.capture, ht.pop, ht.status).aggregate(counts = hl.agg.sum(ht.AC))
    eas_cases = 1716
    eas_controls = 1575
    fin_cases = 381
    fin_controls = 0
    eur_cases = 3188 + 591  + 598
    eur_controls = 5417 + 4017 + 599
    print('get case control counts')
    # change none to 0
    genes = ht.key_by('gene_symbol', 'consequence_category', 'capture', 'pop').select()
    # genes = genes.drop('status')
    genes = genes.distinct() 
    ht = ht.annotate(counts = hl.if_else(hl.is_defined(ht.counts), ht.counts, 0))
    return ht


ht_mac15 = make_rvas_file(ht = ht_in, mac_in = mac_in, mac=15)
ht_mac15.write(OUT, overwrite=True)


