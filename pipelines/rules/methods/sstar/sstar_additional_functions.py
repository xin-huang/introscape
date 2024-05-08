import os

def list_subdirectories(directory):
    subdirectories = []
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path):
            subdirectories.append(item)
    return subdirectories



def process_sstar_1src_output(in_file, out_file):
    """
    Description:
        Helper function for converting output from sstar to BED format files.

    Arguments:
        in_file str: Name of the input file.
        out_file str: Name of the output file.
    """
    df = pd.read_csv(in_file, sep="\t")
    df = df[df['significant'] == True]
    cols = ['chrom', 'start', 'end']
    df.to_csv(out_file, columns=cols, sep="\t", header=False, index=False)


def process_sstar_2src_output(src1_in_file, src2_in_file, threshold_file, src1_out_file, src2_out_file):
    """
    """
    src1_df = pd.read_csv(src1_in_file, sep="\t")
    src2_df = pd.read_csv(src2_in_file, sep="\t")
    threshold_df = pd.read_csv(threshold_file, sep="\t")
    threshold_df = threshold_df[threshold_df['significant'] == True]
    src1_sig_df = pd.merge(src1_df, threshold_df, on=['chrom', 'start', 'end', 'sample'])
    src2_sig_df = pd.merge(src2_df, threshold_df, on=['chrom', 'start', 'end', 'sample'])
    src1_out = open(src1_out_file, 'w')
    src2_out = open(src2_out_file, 'w')
    for i in range(len(src1_sig_df)):
        if (src1_sig_df.iloc[i]['match_pct'] > src2_sig_df.iloc[i]['match_pct']):
            src1_out.write(f'{src1_sig_df.iloc[i]["chrom"]}\t{src1_sig_df.iloc[i]["start"]}\t{src1_sig_df.iloc[i]["end"]}\n')
        if (src1_sig_df.iloc[i]['match_pct'] < src2_sig_df.iloc[i]['match_pct']):
            src2_out.write(f'{src2_sig_df.iloc[i]["chrom"]}\t{src2_sig_df.iloc[i]["start"]}\t{src2_sig_df.iloc[i]["end"]}\n')



def cal_accuracy(true_tracts, inferred_tracts):
    """
    Description:
        Helper function for calculating accuracy.

    Arguments:
        true_tracts str: Name of the BED file containing true introgresssed tracts.
        inferred_tracts str: Name of the BED file containing inferred introgressed tracts.

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """
    truth_tracts = pybedtools.BedTool(true_tracts).sort().merge()
    inferred_tracts =  pybedtools.BedTool(inferred_tracts).sort().merge()

    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

    if float(total_inferred_tracts) == 0: precision = np.nan
    else: precision = true_positives / float(total_inferred_tracts) * 100
    if float(total_true_tracts) == 0: recall = np.nan
    else: recall = true_positives / float(total_true_tracts) * 100

    return precision, recall

