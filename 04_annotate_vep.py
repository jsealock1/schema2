from importlib.metadata import version
from gnomad.utils.vep import process_consequences  

import hail as hl


ht = hl.read_table(HT_VEP)
ht = ht.select('vep')


ht = process_consequences(ht)

# # Case builder function by Konrad, modified by Anne:

print("build anno builder")
PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"] ## same as kyle

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant",
                 'splice_region_variant'] ## kyle only uses msisense 

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"] ## same 

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

# Updated case builder 
# included_loftee_flags should be a list of loftee flags that we will tolerate (not automatically counted as LC)
def annotation_case_builder_updated(worst_csq_for_variant_canonical_expr, 
                                    lof_use_loftee: bool = True, 
                                    included_loftee_flags: list[str] = None,
                                    mis_use_polyphen_and_sift: bool = False, 
                                    mis_use_strict_def: bool = False, 
                                    syn_use_strict_def: bool = False):
    case = hl.case(missing_false=True)
    if lof_use_loftee:
        if included_loftee_flags is None: # If included_loftee_flags is None (WE DON'T WANT TO CONSIDER ANY FLAGS)
            case = (case
                    .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') 
                    .when(worst_csq_for_variant_canonical_expr.lof == 'LC', 'LC'))
        elif len(included_loftee_flags) == 0: # If included_loftee_flags is empty (WE WANT TO FILTER OUT ALL FLAGS)
            case = (case
                    .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') #predicted loss-of-function
                    .when((worst_csq_for_variant_canonical_expr.lof == 'LC') | 
                          (hl.is_defined(worst_csq_for_variant_canonical_expr.lof_flags)), 'LC')) # Low confidence if labelled so, or if there is an additional flag
        elif len(included_loftee_flags) > 0: # If included_loftee_flags has values (WE WANT TO TOLERATE SOME FLAGS)
             case = (case
                    .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') 
                    .when((worst_csq_for_variant_canonical_expr.lof == 'LC') | 
                          (hl.is_defined(worst_csq_for_variant_canonical_expr.lof_flags) & ~hl.set(included_loftee_flags).contains(worst_csq_for_variant_canonical_expr.lof_flags)), 'LC')) # Low confidence if labelled so, or if there is an additional BAD flag
    else:
        case = case.when(hl.set(PLOF_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'pLoF')
    if mis_use_polyphen_and_sift:
        case = (case
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence) &
                      (worst_csq_for_variant_canonical_expr.polyphen_prediction == "probably_damaging") &
                      (worst_csq_for_variant_canonical_expr.sift_prediction == "deleterious"), "damaging_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), "other_missense"))
    else:
        if mis_use_strict_def:
            case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'missense_variant', 'missense')
        else:
            case = case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'missense')
    if syn_use_strict_def:
        case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
    else:
        case = case.when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'synonymous')
    case = case.when(hl.set(OTHER_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'non_coding')
    return case.or_missing()


print('case builder')

# ht = ht.annotate(consequence_category = annotation_case_builder_updated(ht.vep.worst_csq_for_variant_canonical))
possible_flags = ["SINGLE_EXON",
                  "NAGNAG_SITE",
                  "PHYLOCSF_WEAK", 
                  "PHYLOCSF_UNLIKELY_ORF",
                  "NON_CAN_SPLICE"]

these_flags_are_ok = ["NAGNAG_SITE", "NON_CAN_SPLICE","PHYLOCSF_UNLIKELY_ORF", "SINGLE_EXON"]

ht = ht.annotate(consequence_category = annotation_case_builder_updated(ht.vep.worst_csq_for_variant_canonical, 
                                True, 
                                these_flags_are_ok,
                                False, False, False))

print('annotate gerp')
# Annotate in gerp_dist (from vep.worst_csq_for_variant_canonical.lof_info)
# If missing GERP_DIST annotation tag, set score to 0
ht = ht.annotate(gerp_dist = hl.if_else(ht.vep.worst_csq_for_variant_canonical.lof_info.contains('GERP_DIST:'),  # Does it even contain GERP_DIST? 
                                        hl.float64(ht.vep.worst_csq_for_variant_canonical.lof_info.split('GERP_DIST:', 2)[1].split(',', 2)[0]), # If so, extract value and set it as that
                                        hl.float64(0))) # Otherwise, set to 0
# gerp_dist can also be NA if ht.vep.worst_csq_for_variant_canonical.lof_info is NA
ht = ht.annotate(fail_gerp_dist = (hl.is_defined(ht.gerp_dist)) & (ht.gerp_dist < 0)) # Track failures only if there is a GERP_DIST annotation to begin with (ignore NAs)

# Annotate in pass_50_bp_rule (from vep.worst_csq_for_variant_canonical.lof_info)
# If missing, set to false
ht = ht.annotate(pass_50_bp_rule = hl.if_else(ht.vep.worst_csq_for_variant_canonical.lof_info.contains('50_BP_RULE:'),  # Does it even contain 50_BP_RULE?
                                              ht.vep.worst_csq_for_variant_canonical.lof_info.split('50_BP_RULE:', 2)[1].split(',', 2)[0] == "PASS", # If pass, set as true
                                              False)) # Otherwise, set to false
# pass_50_bp_rule can also be NA if ht.vep.worst_csq_for_variant_canonical.lof_info is NA
ht = ht.annotate(fail_50_bp_rule = (hl.is_defined(ht.pass_50_bp_rule)) & (ht.pass_50_bp_rule == False)) # Track failures only if there is a 50_bp_rule annotation to begin with (ignore NAs)

# ht = ht.filter(ht.fail_50_bp_rule & ht.fail_gerp_dist, keep = False)
# create annotation instead of filter
ht = ht.annotate(fail_50_bp_rule_and_gerp_dist = (ht.fail_50_bp_rule & ht.fail_gerp_dist))

# print('Add other annotations')
# ht = ht.annotate(mpc = ht_mpc[ht.locus, ht.alleles])
# ht = ht.annotate(MPC = ht.mpc.MPC)
# ht = ht.drop(ht.mpc)
# ht = ht.annotate(AM = ht_am[ht.locus, ht.alleles].am_pathogenicity)

ht = ht.annotate(
    type = (hl.case()
    .when((hl.len(ht.alleles[0]) == 1) & (hl.len(ht.alleles[1]) == 1), "SNP")
    .when(hl.len(ht.alleles[0]) < hl.len(ht.alleles[1]), "Insertion")
    .when(hl.len(ht.alleles[0]) > hl.len(ht.alleles[1]), "Deletion")
    .or_missing()),
    infrIndel = ((ht.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_insertion") | (ht.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_deletion")),
    )

ht = ht.annotate(gene_symbol = ht.vep.worst_csq_for_variant_canonical.gene_symbol, gene_id = ht.vep.worst_csq_for_variant_canonical.gene_id,
                    transcript = ht.vep.worst_csq_for_variant_canonical.transcript_id)

# csq = ht.aggregate(hl.agg.counter(ht.consequence_category))


ht = ht.annotate(
    isPTV = ht.consequence_category == 'pLoF',
    isLC = ht.consequence_category == 'LC',
    isMIS = ht.consequence_category == 'missense',
    isSYN = ht.consequence_category == 'synonymous',
    isOTH = ht.consequence_category == 'non_coding')

## add in other_splice annotation
ht_os = hl.read_table(HT_OS)
ht = ht.annotate(isOS = hl.is_defined(ht_os[ht.key]))


evaluation_regions = hl.import_locus_intervals(TARGET_INTERVALS)
ht = ht.annotate(eval_reg = hl.is_defined(evaluation_regions[ht.locus]))
ht = ht.filter(ht.isPTV | ht.isLC | ht.isMIS | ht.isSYN | ht.eval_reg, keep = True)




print('write')

ht = ht.repartition(500, shuffle = False)
ht.write(HT_ANNOT, overwrite=True)

