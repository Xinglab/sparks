import pandas as pd
import numpy as np
import os
import argparse
from collections import defaultdict

def import_psi_matrix(input_psi_file, rna_mode=False):

    input_psi_batch = pd.read_csv(input_psi_file, iterator=True, chunksize=10000, index_col=list(range(10)))

    # filter each batch
    input_psi_df_list = []
    mean_filter_range = 0.05
    maxmin_filter_range = 0.1
    for input_psi_block in input_psi_batch:
        # filter by psi mean
        psi_mean = input_psi_block.mean(axis=1)
        mean_filtered_block = input_psi_block[(psi_mean > mean_filter_range) & (psi_mean < 1 - mean_filter_range)]

        # filter by max min
        max_min = mean_filtered_block.max(axis=1) - mean_filtered_block.min(axis=1)
        maxmin_filtered_block = mean_filtered_block[max_min > maxmin_filter_range]

        # add it for later
        input_psi_df_list.append(maxmin_filtered_block)
    input_psi_df = pd.concat(input_psi_df_list)


    if not rna_mode:
        # update column names - the sample names are being used with the last run identifier (ABCD)
        # so this need to be trimmed in order to be matched with the MAF file
        sample_list = input_psi_df.columns.tolist()
        sample_list_trimmed = [sample[:-1] for sample in sample_list]
        input_psi_df.columns = sample_list_trimmed

    return input_psi_df

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(description="")

    input_args = parser.add_argument_group("Input Arguments")
    input_args.add_argument("input_PSI_file", help="Calculated PSI matrix")
    input_args.add_argument("output_dir", help="output dir to store select MAF file")

    input_args.add_argument('--rna', help="RNA mode where sample matching is not necessary",
                            action="store_true")
    input_args.add_argument("--cancer_type", help="cancer type - e.g. TCGA-LAML - default = TCGA-BRCA",
                            action="store", default="TCGA-BRCA")

    args = parser.parse_args()

    input_psi_file = args.input_PSI_file
    output_dir = args.output_dir
    rna_mode = args.rna
    study = args.cancer_type

    splice_type = input_psi_file.split('/')[-1].split('.')[1]

    print("Importing PSI matrix")
    input_psi_df = import_psi_matrix(input_psi_file, rna_mode)

    # write the filtered PSI df
    print("Writing PSI matrix")
    filtered_psi_file = '%s/filtered_psi.matrix.%s.%s.csv' % (output_dir, study, splice_type)
    input_psi_df.to_csv(filtered_psi_file)

    print("PSI file stored in %s" % filtered_psi_file)