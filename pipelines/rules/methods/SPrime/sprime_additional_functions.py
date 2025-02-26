#import pybedtools
import numpy as np
import pandas as pd

def create_map_file(x, y, map_file = "sim.map"):
    with open(map_file, 'w') as f:
        f.write("1\t.\t0.0\t0\n")
        f.write(f"1\t.\t{x}\t{y}\n")

#-----------------------------------------------------------------------------------------------------------------------
# process_sprime_output fn from https://github.com/admixVIE/sstar-analysis/blob/main/utils/utils.py
#-----------------------------------------------------------------------------------------------------------------------

def process_sprime_output(in_file, out_file):
    """
    Description:
        Helper function for converting output from SPrime to BED format.

    Arguments:
        in_file str: Name of the input file.
        out_file str: Name of the output file.
    """
    # read in the data - the SPrime output
    df = pd.read_csv(in_file, delimiter="\t")

    # drop columns that are not needed
    df2 = df.drop(['ID', 'REF', 'ALT', 'ALLELE'], axis=1)

    # add columns START and END with the highest ans lowest position of each chromosome, segment and score
    df2['START'] = df2.groupby(['CHROM', 'SCORE', 'SEGMENT'])['POS'].transform(min)
    df2['END'] = df2.groupby(['CHROM', 'SCORE', 'SEGMENT'])['POS'].transform(max)

    # group by chromosome, segment and score - drop the column position
    df3 = df2.loc[df2.groupby(["CHROM", "SCORE", "SEGMENT"])["POS"].idxmax()]
    df4 = df3.drop(['POS'], axis=1)

    # get the right order (for the bed file)
    df_final = df4[['CHROM','START','END','SEGMENT','SCORE']].sort_values(by=['START', 'SEGMENT'])

    np.savetxt(out_file, df_final.values, fmt='%s', delimiter='\t')


#-----------------------------------------------------------------------------------------------------------------------
# rm cal_accuracy fn here, & call instead from pipelines/rules/commons/evaluate_utils.py
#-----------------------------------------------------------------------------------------------------------------------

#def cal_accuracy(true_tracts, inferred_tracts):
#    """
#    Description:
#        Helper function for calculating accuracy.
#
#    Arguments:
#        true_tracts str: Name of the BED file containing true introgresssed tracts.
#        inferred_tracts str: Name of the BED file containing inferred introgressed tracts.
#
#    Returns:
#        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
#        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
#    """
#    truth_tracts = pybedtools.BedTool(true_tracts).sort().merge()
#    inferred_tracts =  pybedtools.BedTool(inferred_tracts).sort().merge()
#
#    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
#    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
#    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])
#
#    if float(total_inferred_tracts) == 0: precision = np.nan
#    else: precision = true_positives / float(total_inferred_tracts) * 100
#    if float(total_true_tracts) == 0: recall = np.nan
#    else: recall = true_positives / float(total_true_tracts) * 100
#
#    return precision, recall

