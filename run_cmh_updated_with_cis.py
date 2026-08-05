import hail as hl
import argparse
import math

# ---------------------
# Hail init
# ---------------------

hl.init(
    default_reference='GRCh38'
)

def compute_CMH_OR_with_CI(case_carriers, case_non_carriers,
                            control_carriers, control_non_carriers,
                            alpha=0.05):
    """
    Returns (OR, SE_log_OR, CI_lower, CI_upper) using the
    Robins-Breslow-Greenland variance estimator for the CMH OR.
    """

    def P_term(a, b, c, d, t):
        # P_i = (a_i + d_i) / t_i
        return (a + d) / t

    def Q_term(a, b, c, d, t):
        # Q_i = (b_i + c_i) / t_i
        return (b + c) / t

    def R_term(a, d, t):
        # R_i = (a_i * d_i) / t_i
        return (a * d) / t

    def S_term(b, c, t):
        # S_i = (b_i * c_i) / t_i
        return (b * c) / t

    n1 = hl.zip(case_carriers, case_non_carriers).map(lambda ab: ab[0] + ab[1])
    n2 = hl.zip(control_carriers, control_non_carriers).map(lambda cd: cd[0] + cd[1])
    t  = hl.zip(n1, n2).map(lambda nn: nn[0] + nn[1])

    # Precompute per-stratum R_i and S_i (reused for OR and variance)
    R = hl.zip(case_carriers, control_non_carriers, t).map(
        lambda tup: R_term(*tup)
    )
    S = hl.zip(case_non_carriers, control_carriers, t).map(
        lambda tup: S_term(*tup)
    )

    sum_R = hl.sum(R)
    sum_S = hl.sum(S)
    OR = sum_R / sum_S

    # --- Robins-Breslow-Greenland variance of log(OR) ---
    # Var(log OR) = sum(P_i * R_i) / (2 * sum(R_i)^2)
    #             + sum(P_i * S_i + Q_i * R_i) / (2 * sum(R_i) * sum(S_i))
    #             + sum(Q_i * S_i) / (2 * sum(S_i)^2)

    P = hl.zip(case_carriers, case_non_carriers,
               control_carriers, control_non_carriers, t).map(
        lambda tup: P_term(*tup)
    )
    Q = hl.zip(case_carriers, case_non_carriers,
               control_carriers, control_non_carriers, t).map(
        lambda tup: Q_term(*tup)
    )

    PR = hl.sum(hl.zip(P, R).map(lambda pr: pr[0] * pr[1]))
    QS = hl.sum(hl.zip(Q, S).map(lambda qs: qs[0] * qs[1]))
    PS_QR = hl.sum(
        hl.zip(P, S, Q, R).map(lambda tup: tup[0] * tup[1] + tup[2] * tup[3])
    )

    var_log_OR = (
          PR       / (2 * sum_R ** 2)
        + PS_QR    / (2 * sum_R * sum_S)
        + QS       / (2 * sum_S ** 2)
    )

    SE_log_OR = hl.sqrt(var_log_OR)

    # z* for two-sided CI (default alpha=0.05 → z=1.96)
    z = -math.log(alpha / 2)          # ≈ 1.96 for alpha=0.05; pure Python scalar
    # or hardcode: z = 1.959964

    CI_lower = OR * hl.exp(-z * SE_log_OR)
    CI_upper = OR * hl.exp( z * SE_log_OR)

    return OR, SE_log_OR, CI_lower, CI_upper




def main(args):
    print('read mt')

    mt = hl.read_matrix_table(args.mt)


    print('create counts mt from count data')

    count_tbl = hl.read_table(args.counts)

    # Harmonize CASE/CTRL
    count_tbl = count_tbl.annotate(
        case_con = hl.if_else(count_tbl.status == 'schizophrenia', 'CASE', 'CTRL')
    )

    # Harmonize populations: EUR→NFE, FIN→FIN, EAS→EAS
    count_tbl = count_tbl.annotate(
        pop2 = hl.if_else(count_tbl.pop == 'EUR', 'NFE', count_tbl.pop)
    )
    count_tbl = count_tbl.annotate(
        pop3 = hl.if_else(count_tbl.pop2 == 'FIN', 'FIN', count_tbl.pop2)
    )
    count_tbl = count_tbl.annotate(
        pop4 = hl.if_else(count_tbl.pop3 == 'EAS', 'EAS', count_tbl.pop3)
    )

    # group2: pop_capture_caseCon, e.g. "EAS_Agilent_CASE"
    count_tbl = count_tbl.annotate(
        group2 = count_tbl.pop4 + "_" + count_tbl.capture + "_" + count_tbl.case_con
    )

    count_tbl = count_tbl.key_by()
    count_tbl = count_tbl.select(
        'group2', 'case_con', 'gene_symbol', 'consequence_category', 'counts'
    )

    # Matrix with counts as entries
    count_mt = count_tbl.to_matrix_table(
        row_key=['gene_symbol', 'consequence_category'],
        col_key=['group2', 'case_con']
    )

    # Collapse columns by group2 to get n_carriers_add per gene × group
    count_mt = count_mt.group_cols_by(count_mt.group2).aggregate(
        n_carriers_add = hl.agg.sum(count_mt.counts)
    )
    count_mt = count_mt.annotate_cols(
        case_con = count_mt.group2.split("_")[-1]
    )

    # Optional consequence filters
    if args.ptv_and_mis:
        print('filtering counts to ptv and missense')
        count_mt = count_mt.filter_rows(
            (count_mt.consequence_category == 'ptv') |
            (count_mt.consequence_category == 'mis_mean_rank93')
        )

    if args.ptv_only:
        print('filtering counts to ptv only')
        count_mt = count_mt.filter_rows(count_mt.consequence_category == 'ptv')

    if args.mis_only:
        print('filtering counts to misrank93 only')
        count_mt = count_mt.filter_rows(count_mt.consequence_category == 'mis_mean_rank93')
    
    if args.ptv_and_mpc3:
        print('filtering counts to ptv and missense mpc >3')
        count_mt = count_mt.filter_rows(
            (count_mt.consequence_category == 'ptv') |
            (count_mt.consequence_category == 'mpc3')
        )
    
    if args.ptv_and_mpc2:
        print('filtering counts to ptv and missense mpc 2-3')
        count_mt = count_mt.filter_rows(
            (count_mt.consequence_category == 'ptv') |
            (count_mt.consequence_category == 'mpc2_3')
        )

    if args.mpc3:
        print('filtering counts to missense mpc >3')
        count_mt = count_mt.filter_rows(
            (count_mt.consequence_category == 'mpc3')
        )

    if args.synonymous:
        print('filtering counts to synonymous')
        count_mt = count_mt.filter_rows(
            (count_mt.consequence_category == 'synonymous')
        )



    count_mt = count_mt.key_rows_by('gene_symbol')
    count_mt = count_mt.drop('consequence_category')


    # ---------------------
    # External sample totals (Agilent etc., including EAS)
    # These are per-group2 totals (not per-gene)
    # ---------------------

    counts_add = hl.literal({
        'EAS_Agilent_CASE': 1716,
        'EAS_Agilent_CTRL': 1575,
        'FIN_Agilent_CASE': 381,
        'FIN_Agilent_CTRL': 0,
        'NFE_Agilent_CASE': 4377,
        'NFE_Agilent_CTRL': 10033
    })


    print('Build per-gene EAS Agilent counts from count_mt')

    # Keep only EAS Agilent CASE/CTRL columns
    eas_counts = count_mt.filter_cols(
        (count_mt.group2 == 'EAS_Agilent_CASE') |
        (count_mt.group2 == 'EAS_Agilent_CTRL')
    )

    # Just need gene + group2 + n_carriers_add at entry level
    eas_entries = eas_counts.select_entries('n_carriers_add').entries()

    # Aggregate to per-gene CASE/CTRL carrier counts
    eas_tbl = eas_entries.group_by(eas_entries.gene_symbol).aggregate(
        eas_case_carriers = hl.agg.sum(
            hl.if_else(eas_entries.group2 == 'EAS_Agilent_CASE',
                    eas_entries.n_carriers_add, 0)
        ),
        eas_ctrl_carriers = hl.agg.sum(
            hl.if_else(eas_entries.group2 == 'EAS_Agilent_CTRL',
                    eas_entries.n_carriers_add, 0)
        )
    ).key_by('gene_symbol')


    print('filter mt from manifest')

    manifest_ht = hl.import_table(
        args.manifest,
        delimiter="\t",
        impute=True,
        key="s"
    )

    # Keep only samples that appear in manifest and are not "exclude"
    mt = mt.filter_cols(
        hl.if_else(
            hl.is_defined(manifest_ht[mt.s]),
            manifest_ht[mt.s]['status'] != "exclude",
            False
        ),
        keep=True
    )

    # Annotate columns: status → CASE/CTRL; pop; capture; group keys
    mt = mt.annotate_cols(
        case_con_raw = manifest_ht[mt.s]['status'],
        case_con     = hl.if_else(manifest_ht[mt.s]['status'] == 'schizophrenia',
                                "CASE", "CTRL"),
        pop          = manifest_ht[mt.s]['pop'],
        chip         = manifest_ht[mt.s]['capture']
    )

    mt = mt.annotate_cols(
        group  = mt.pop + "_" + mt.chip,
        group2 = mt.pop + "_" + mt.chip + "_" + mt.case_con
    )

    # Basic group bookkeeping (optional; debug)
    groups = mt.aggregate_cols(hl.agg.collect_as_set(mt.group))
    group2_counts = mt.aggregate_cols(hl.agg.counter(mt.group2))

    for g in groups:
        if g + "_CASE" not in group2_counts:
            group2_counts[g + "_CASE"] = 0
        if g + "_CTRL" not in group2_counts:
            group2_counts[g + "_CTRL"] = 0

    group2_counts = hl.literal(group2_counts)

    # Repartition / persist for speed
    if mt.n_partitions() > 200:
        mt = mt.repartition(200, shuffle=False)

    mt = mt.persist()


    print('Aggregate internal carriers within MT by gene × group2')


    m = mt  # alias

    res = m.group_cols_by(m.group2).aggregate(
        n_carriers = hl.agg.sum(m.agg > 0),
        n_total    = hl.agg.count()
    ).repartition(200, shuffle=False).persist()

    # Internal non-carriers
    res = res.annotate_entries(
        n_non_carriers = res.n_total - res.n_carriers
    )

    # Recover CASE/CTRL from group2
    res = res.annotate_cols(
        case_con = res.group2.split("_")[-1]
    )


    # ---------------------
    print('Add external counts (FIN/NFE/EAS, etc.) from count_mt and counts_add')
    # ---------------------

    # Add external carriers (if present in count_mt)
    res = res.annotate_entries(
        n_carriers_add_ext = hl.or_else(
            count_mt[res.row_key, res.col_key].n_carriers_add,
            0
        )
    )

    # Total external samples per group2 (from counts_add)
    res = res.annotate_cols(
        n_total_add_ext = hl.or_else(counts_add.get(res.group2), 0)
    )

    # Combined totals per cell: internal + external
    res = res.transmute_entries(
        n_carriers_tot     = hl.int64(res.n_carriers) +
                            hl.int64(res.n_carriers_add_ext),
        n_non_carriers_tot = hl.int64(res.n_non_carriers) +
                            (hl.int64(res.n_total_add_ext) -
                            hl.int64(res.n_carriers_add_ext))
    )


    # ---------------------
    print('Build CMH strata from combined totals')
    # base = everything before _CASE/_CTRL in group2
    # ---------------------

    res = res.annotate_cols(
        base    = hl.delimit(res.group2.split("_")[:-1], "_"),
        is_case = (res.case_con == "CASE")
    )

    # Collect per-row arrays of structs for CASE / CTRL (internal + non-EAS external)
    res_counts = res.annotate_rows(
        case_arr = hl.agg.collect(
            hl.struct(
                base       = res.base,
                n_carriers = hl.if_else(res.is_case,  res.n_carriers_tot,     0),
                n_noncar   = hl.if_else(res.is_case,  res.n_non_carriers_tot, 0)
            )
        ),
        ctrl_arr = hl.agg.collect(
            hl.struct(
                base       = res.base,
                n_carriers = hl.if_else(~res.is_case, res.n_carriers_tot,     0),
                n_noncar   = hl.if_else(~res.is_case, res.n_non_carriers_tot, 0)
            )
        )
    ).rows()


    print('add EAS Agilent as an extra CMH stratum per gene')

    # Attach per-gene EAS counts (may be missing for many genes)
    res_counts = res_counts.annotate(
        eas_info = eas_tbl[res_counts.gene_symbol]
    )

    # 1) carriers: treat missing eas_info as 0 carriers
    res_counts = res_counts.annotate(
        eas_case_carriers = hl.int64(
            hl.or_else(res_counts.eas_info.eas_case_carriers, 0)
        ),
        eas_ctrl_carriers = hl.int64(
            hl.or_else(res_counts.eas_info.eas_ctrl_carriers, 0)
        )
    )

    # 2) noncarriers: always full totals – carriers subtracted
    res_counts = res_counts.annotate(
        eas_case_noncarriers = hl.int64(counts_add['EAS_Agilent_CASE']) -
                            res_counts.eas_case_carriers,
        eas_ctrl_noncarriers = hl.int64(counts_add['EAS_Agilent_CTRL']) -
                            res_counts.eas_ctrl_carriers
    )

    # 3) Append EAS_Agilent as a struct into case_arr / ctrl_arr for EVERY gene
    res_counts = res_counts.annotate(
        case_arr = res_counts.case_arr.append(
            hl.struct(
                base       = 'EAS_Agilent',
                n_carriers = res_counts.eas_case_carriers,
                n_noncar   = res_counts.eas_case_noncarriers
            )
        ),
        ctrl_arr = res_counts.ctrl_arr.append(
            hl.struct(
                base       = 'EAS_Agilent',
                n_carriers = res_counts.eas_ctrl_carriers,
                n_noncar   = res_counts.eas_ctrl_noncarriers
            )
        )
    )

    # ---------------------
    print('Summarize by base → CMH input arrays')
    # ---------------------

    def _sum_by_base(arr, field):
        return hl.group_by(lambda x: x.base, arr).map_values(
            lambda xs: hl.sum(xs.map(lambda y: y[field]))
        )

    res_counts = res_counts.annotate(
        case_car_dict    = _sum_by_base(res_counts.case_arr,  'n_carriers'),
        case_noncar_dict = _sum_by_base(res_counts.case_arr,  'n_noncar'),
        ctrl_car_dict    = _sum_by_base(res_counts.ctrl_arr,  'n_carriers'),
        ctrl_noncar_dict = _sum_by_base(res_counts.ctrl_arr,  'n_noncar'),
        bases = hl.sorted(
            hl.set(res_counts.case_arr.map(lambda x: x.base))
            .union(hl.set(res_counts.ctrl_arr.map(lambda x: x.base)))
        )
    )

    # Turn dicts into aligned arrays for CMH (EAS_Agilent included via case_arr/ctrl_arr)
    res_counts = res_counts.annotate(
        case_carriers        = res_counts.bases.map(
            lambda b: hl.or_else(res_counts.case_car_dict.get(b), 0)
        ),
        case_non_carriers    = res_counts.bases.map(
            lambda b: hl.or_else(res_counts.case_noncar_dict.get(b), 0)
        ),
        control_carriers     = res_counts.bases.map(
            lambda b: hl.or_else(res_counts.ctrl_car_dict.get(b), 0)
        ),
        control_non_carriers = res_counts.bases.map(
            lambda b: hl.or_else(res_counts.ctrl_noncar_dict.get(b), 0)
        )
    )


    # Totals per gene across all strata (including EAS_Agilent)
    res_counts = res_counts.annotate(
        total_case_carriers        = hl.sum(res_counts.case_carriers),
        total_case_non_carriers    = hl.sum(res_counts.case_non_carriers),
        total_control_carriers     = hl.sum(res_counts.control_carriers),
        total_control_non_carriers = hl.sum(res_counts.control_non_carriers)
    )


    # ---------------------
    print('Run CMH per gene')
    # ---------------------


    res_counts = res_counts.drop(
        'eas_info',
        'eas_case_carriers',
        'eas_ctrl_carriers',
        'eas_case_noncarriers',
        'eas_ctrl_noncarriers'
    )

    cmh_or, cmh_se, cmh_ci_lower, cmh_ci_upper = (
        compute_CMH_OR_with_CI(
            res_counts.case_carriers,
            res_counts.case_non_carriers,
            res_counts.control_carriers,
            res_counts.control_non_carriers
        )
    )

    res_counts = res_counts.annotate(
        RES        = hl.cochran_mantel_haenszel_test(
                        res_counts.case_carriers,
                        res_counts.case_non_carriers,
                        res_counts.control_carriers,
                        res_counts.control_non_carriers
                    ),
        OR         = cmh_or,
        SE_log_OR  = cmh_se,
        CI_lower   = cmh_ci_lower,
        CI_upper   = cmh_ci_upper
    ).flatten()

    # ---------------------
    print('Export')
    # ---------------------

    res_counts.export(args.out)



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mt",
        help = "Path to starting VEP-counts MT (computed in 02_VEP-counts-export.py)",
        type = str,
        required = True
    )
    parser.add_argument(
        "--manifest",
        help = "Path to tab-delimited manifest keyed by samples ('s') with relevant annotations: 'FILTER', 'CASECON', 'POP', 'CHIP'",
        type = str,
        required = True
    )
    parser.add_argument(
        "--counts",
        help = "Path to additional counts hail table",
        type = str,
        required = False
    )
    parser.add_argument(
        "--out",
        help = "Path to output directory",
        type = str,
        required = True
    )
    parser.add_argument(
        "--ptv_only",
        help = "Use only PTVs",
        action = 'store_true'
    )  
    parser.add_argument(
        "--ptv_and_mis",
        help = "Use PTVs and Missense variants",
        action = 'store_true'
    )  
    parser.add_argument(
        "--mis_only",
        help = "Use Missense variants",
        action = 'store_true'
    )  

    parser.add_argument(
        "--ptv_and_mpc3",
        help = "Use PTVs and Missense variants",
        action = 'store_true'
    )  

    parser.add_argument(
        "--ptv_and_mpc2",
        help = "Use PTVs and Missense variants",
        action = 'store_true'
    )  

    parser.add_argument(
        "--mpc3",
        help = "Use PTVs and Missense variants",
        action = 'store_true'
    )  
    parser.add_argument(
        "--synonymous",
        help = "Use syn variants",
        action = 'store_true'
    )  
    args = parser.parse_args()

    main(args)
