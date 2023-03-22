import pandas as pd
import argparse
import os
import numpy as np

def process_mats(input_file,
                 event_type,
                 output_dir, 
                 study,
                 one_group_mode = False):
    # import mats file
    input_mats = pd.read_csv(input_file, sep='\t', index_col=0)

    clean_mats_list = []
    for idx, event_entry in input_mats.iterrows():

        # concatenate relevant information to make ID
        event_id = ":".join([str(x) for x in event_entry[0:10]] + [event_type])

        # calculate pulled psi values
        # pull counts together
        inc_count_group1 = sum([int(x) for x in event_entry['IJC_SAMPLE_1'].split(',')])
        skip_count_group1 = sum([int(x) for x in event_entry['SJC_SAMPLE_1'].split(',')])
        if not one_group_mode:
            inc_count_group2 = sum([int(x) for x in event_entry['IJC_SAMPLE_2'].split(',')])
            skip_count_group2 = sum([int(x) for x in event_entry['SJC_SAMPLE_2'].split(',')])

        # query length for psi calculation
        inc_len = event_entry['IncFormLen']
        skip_len = event_entry['SkipFormLen']

        try:
            # calculate new psi
            psi_group1 = (inc_count_group1 / inc_len) / ((inc_count_group1 / inc_len) + (skip_count_group1 / skip_len))
            if not one_group_mode:
                psi_group2 = (inc_count_group2 / inc_len) / ((inc_count_group2 / inc_len) + (skip_count_group2 / skip_len))

                delta_psi_pulled = psi_group1 - psi_group2
            else:  # if one group mode - do not calculate this 
                psi_group2 = np.nan
                delta_psi_pulled = np.nan
        except: 
            continue

        # calculate average count
        if not one_group_mode:
            inc_count = [int(x) for x in event_entry['IJC_SAMPLE_1'].split(',') + event_entry['IJC_SAMPLE_2'].split(',')]
            skip_count = [int(x) for x in event_entry['SJC_SAMPLE_1'].split(',') + event_entry['SJC_SAMPLE_2'].split(',')]
        else: 
            inc_count = [int(x) for x in event_entry['IJC_SAMPLE_1'].split(',')]
            skip_count = [int(x) for x in event_entry['SJC_SAMPLE_1'].split(',')]

        total_count = [x + y for x, y in zip(inc_count, skip_count)]
        avg_count = sum(total_count) / len(total_count)

        # merge raw count values for downstream diagnostics count analysis
        count_values = ",".join([str(x) for x in total_count])

        # merge raw psi values for downstream processing
        if not one_group_mode:
            psi_values = ",".join([event_entry['IncLevel1'], event_entry['IncLevel2']])
        else:
            psi_values = event_entry['IncLevel1']

        # extract relevant information
        beta = event_entry['IncLevelDifference']
        pval = event_entry['PValue']
        fdr = event_entry['FDR']


        clean_mats_entry = pd.DataFrame(data=[event_id, beta, pval, fdr, avg_count, psi_values, count_values,
                                              inc_len, skip_len, psi_group1, psi_group2, delta_psi_pulled])
        clean_mats_list.append(clean_mats_entry)

    # append to make a dataframe
    clean_mats_df = pd.concat(clean_mats_list, axis=1).transpose()
    clean_mats_df.columns = ['event', 'beta', 'pval', 'fdr', 'avg_count', 'psi_values', 'count_values',
                             'inc_len', 'skip_len', 'pulled_psi_1', 'pulled_psi_2', 'pulled_delta_psi']
    # save the output file
    output_file = "%s/%s.%s.MATS_df.txt" % (output_dir, study, event_type)
    clean_mats_df.to_csv(output_file, sep='\t', index=False)

    return(output_file)

if __name__ == "__main__":

    script_dir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(description="Cleanup MATS results for downstream processing")

    input_args = parser.add_argument_group("Input Arguments")
    input_args.add_argument("input_dir", help="directory where the MATS result are stored")
    input_args.add_argument("output_dir", help="output dir to store filtered MATS df file")
    input_args.add_argument("study", help="Study to annotate the output file")
    input_args.add_argument('--one_group', help="One group mode for rMATS run without stats mode", action="store_true")
    input_args.add_argument('--SE', help="Run this script for SE only", action="store_true")

    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    study = args.study
    one_group_mode = args.one_group
    se_mode = args.SE

    # make analysis dir if it does not exist
    if not os.path.exists(output_dir):  
        os.makedirs(output_dir)

    if se_mode:
        event_type_list = ['SE']
    else: 
        event_type_list = ['SE', 'A3SS', 'A5SS', 'RI']

    # process 
    for event_type in event_type_list:
        print("Processing MATS data for %s - %s" % (study, event_type))
        
        # process mats file 
        input_mats_file = "%s/%s.MATS.JC.txt" % (input_dir, event_type)

        output_file = process_mats(input_mats_file, event_type, output_dir, study, one_group_mode)

        print("Output file stored in %s" % output_file)








