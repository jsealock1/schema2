"""
Get counts for specific VEP annotations and write 
"""

print('import packages')
import hail as hl
import logging
import argparse
from datetime import date

print('logging')
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def main(args):
    # Initialize
    print('init')
    hl.init(tmp_dir = args.tmp)
    hl.default_reference('GRCh38')
        
    # Read in MT
    logger.info(f"Reading in MT ({args.mt})...")
    mt = hl.read_matrix_table(args.mt)

    logger.info(f"Starting (variants, samples): {mt.count()}")

    # Read in VEP HT
    ht = hl.read_table(args.vep_ht)

    # Filter to variants in VEP HT (usually already filtered to consequence category in pLoF, other_missense, damaging_missense, synonymous)
    mt = mt.filter_rows(hl.is_defined(ht[mt.locus, mt.alleles]), keep = True)

    # Adjust consequence annotations here - provided by Julia
    ht = ht.filter(ht.fail_50_bp_rule_and_gerp_dist, keep=False)
    ht = ht.annotate(consequence = hl.if_else((ht.consequence_category == "pLoF") | (ht.isOS), 'ptv',
                                hl.if_else(ht.consequence_category == "LC", 'ptv_lc',
                                hl.if_else((ht.consequence_category == "synonymous") & (ht.isOS==False), 'synonymous',
                                hl.if_else(((ht.mis_mean_rank >= 93) & (ht.isOS==False)),  'mis_mean_rank93', 
                                'none')))))


    # Annotate in relevant fields from VEP HT
    mt = mt.annotate_rows(gene_symbol = ht[mt.locus, mt.alleles].gene_symbol,
                            consequence = ht[mt.locus, mt.alleles].consequence)

    f = [args.file_prefix] # Running identifier/run details

    ## Next, consequence_category
    if args.cons_cat:
        cons_cat = args.cons_cat.replace(" ", "").split(",") # Parse command-line input
        logger.info(f"Filtering to {cons_cat}...")
        mt = mt.filter_rows(hl.set(cons_cat).contains(mt.consequence), keep = True)
        f.append(f"{'-'.join(cons_cat)}")

    ## Next, AC threshold
    if args.ac == 1:
        logger.info(f"Filtering to singletons...")
        if args.ac_file:
            ac = hl.read_table(args.ac_file)
            ac = ac.filter(ac.MAC == 1, keep=True)
            mt = mt.filter_rows(hl.is_defined(ac[mt.locus, mt.alleles]))
        else:
            mt = mt.filter_rows(mt.variant_qc.AC[1] == 1, keep = True)
        f.append(f"singletons")
    elif args.ac:
        logger.info(f"Filtering to AC <= {args.ac}...")
        if args.ac_file:
            ac = hl.read_table(args.ac_file)
            ac = ac.filter(ac.MAC <= args.ac, keep=True)
            mt = mt.filter_rows(hl.is_defined(ac[mt.locus, mt.alleles]))
        else:
            mt = mt.filter_rows(mt.variant_qc.AC[1] <= args.ac, keep = True)
        f.append(f"AC{args.ac}")
    

    logger.info(f"Generating counts...")
    mt_agg = mt.group_rows_by(mt.gene_symbol).aggregate(
        agg = hl.agg.count_where(mt.GT.is_non_ref())
    )

    logger.info(f"Writing to {args.out + '_'.join(f) + '_counts.mt'}...")
    print('1')
    mt_agg.repartition(200, shuffle = False).write(args.out + '_'.join(f) + '_counts.mt', overwrite = True)
    print('2')

    logger.info(f"Copying log to {args.out + 'logs/' + str(date.today()) + '_' + '_'.join(f) + '.log'}...")
    hl.copy_log(f"{args.out + 'logs/' + str(date.today()) + '_' + '_'.join(f) + '.log'}")




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mt",
        help = "Path to input MT (computed in 00_filter-to-coding.py)",
        type = str,
        required = True
    )
    parser.add_argument(
        "--vep_ht",
        help = "Path to VEP-annotation HT for input MT (computed in 01_vep-annotate_[run].py)",
        type = str,
        required = True
    )
    parser.add_argument(
        "--gnomAD_AC_thresh",
        help = "gnomAD AC threshold to filter to (keep entries <= thresh)",
        type = int,
        required = False
    )   
    parser.add_argument(
        "--RGC_AC_thresh",
        help = "RGC AC threshold to filter to (keep entries <= thresh)",
        type = int,
        required = False
    )  
    parser.add_argument(
        "--cons_cat",
        help = "Comma-separated consequence categories to filter to: pLoF, other_missense, damaging_missense, missense, synonymous",
        type = str,
        required = False
    )
    parser.add_argument(
        "--ac",
        help = "AC threshold to filter to (remove entries < thresh), usually 5 or 10",
        type = int,
        required = False
    )
    parser.add_argument(
        "--ac_file",
        help = "Path to combined allele counts (including iPSYCH/trios/UK10K)",
        type = str,
        required = False
    )
    parser.add_argument(
        "--mpc",
        help = "MPC threshold to filter to (remove entries < thresh), usually 2 (with AM 0.98)",
        type = float,
        required = False
    )
    parser.add_argument(
        "--am",
        help = "AlphaMissense threshold to filter to (>= thresh), usually 0.98 (with MPC 2)",
        type = float,
        required = False
    )
    parser.add_argument(
        "--misfitS",
        help = "misfitS threshold to filter to (> thresh), usually 0.03",
        type = float,
        required = False
    )      
    parser.add_argument(
        "--union_missense",
        help = "Flag for whether missense filters should be considered as union (intersection by default) ",
        action = 'store_true'
    )   
    parser.add_argument(
        "--non_gnomAD_psych",
        help = "Flag for filtering to non-gnomAD-psych variants",
        action = 'store_true'
    )
    parser.add_argument(
        "--out",
        help = "Path to output directory",
        type = str,
        required = True
    )
    parser.add_argument(
        "--file_prefix",
        help = "File output prefix",
        type = str,
        required = True
    )
    parser.add_argument(
        "--tmp",
        help = "Path to temp bucket",
        type = str,
        required = True
    )

    args = parser.parse_args()

    main(args)