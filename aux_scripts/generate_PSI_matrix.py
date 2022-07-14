import sys
import os
import pandas as pd
import numpy as np
import argparse
import time
import math
import pandas as pd

def calculate_PSI(inc_list, skip_list, inclusion_exon_size, skipping_exon_size, event):
    psi_list = []
    for i in range(len(inc_list)):
        inclusion_value = float(inc_list[i])
        skipping_value = float(skip_list[i])
        try:
            if inclusion_value + skipping_value >= 10:
                psi_value = round((inclusion_value / inclusion_exon_size) / (
                            (inclusion_value / inclusion_exon_size) + (skipping_value / skipping_exon_size)), 4)
            else:
                psi_value = np.nan
        except ZeroDivisionError:
            psi_value = np.nan

        psi_list.append(psi_value)
    return psi_list


def read_sample(sample_file):
    sample_list = []
    csv_indicator = 0
    with open(sample_file, 'r') as f:
        for line in f:  # detect if csv
            if len(line.split(',')) > 1:  # if csv
                samples = line.rstrip().split(',')
                for sample in samples:
                    sample_name = sample.split('/')[-1].split('.bam')[0]
                    sample_list.append(sample_name)
                csv_indicator = 1

            if csv_indicator == 0:  # if new line
                sample = line.rstrip()
                sample_name = sample.split('/')[-1].split('.bam')[0]
                sample_list.append(sample_name)
    return sample_list


def generate_psi_matrix(from_gtf, jc_file, sample_file, output_dir, tissue, event_type,
                        two_group_mode=False, second_sample_file="",
                        coverage_filter=False, coverage_cutoff=10,
                        NA_filter=False, NA_cutoff=0.9,
                        split_by_chromosome=False,
                        chromosome=""):
    line_count = 0
    df = pd.DataFrame()
    new_vector_list = []
    if not two_group_mode:
        iterator = [0]
    elif two_group_mode:
        iterator = [0, 2]

    num_below_coverage = 0
    num_below_NA = 0
    # make PSI matrix without SNP information

    with open(jc_file, 'r') as f:
        for line in f:

            if line_count == 0:
                line_count += 1
                continue

            line_count += 1
            if line_count % 1000 == 0:
                print(("processed %i event " % line_count))
                print(("number of events = %i" % len(inc)))

            inc = []
            skip = []
            for i in iterator:  # iterate twice for two sample mode
                elements = line.split('\t')
                event_ID = str(elements[0])
                inc += elements[i + 1].split(',')
                skip += elements[i + 2].split(',')

            if coverage_filter:  # if coverage filter is enabled
                coverage = [int(x) + int(y) for x, y in zip(inc, skip)]  # if covered at all
                # print coverage
                if np.mean(coverage) < coverage_cutoff:  # if below coverage cutoff
                    num_below_coverage += 1
                    continue  # skip this event

            inc_exon_size = float(elements[5])
            skip_exon_size = float(elements[6])
            psi_list = calculate_PSI(inc, skip, inc_exon_size, skip_exon_size, event_type)

            if NA_filter:  # if na filter is enabled
                na_count = sum(np.isnan(psi_list))
                nan_cutoff_number = math.ceil(NA_cutoff * len(psi_list))
                if na_count > nan_cutoff_number:  # if more samples than cutoff are missing
                    num_below_NA += 1
                    continue

            # if line_count % 100 == 0:
            # 	print inc
            # 	print skip
            # 	print psi_list

            psi_df = pd.DataFrame({event_ID: psi_list})  # generate PSI vector
            new_vector_list.append(psi_df)

        # print psi_df
    print((line_count, len(new_vector_list), num_below_coverage, num_below_NA))

    df = pd.concat(new_vector_list, axis=1)  # merge the matrix and vector

    # get sample
    # TODO - for now in two sample mode, it requires the sample files to be concatenated - fix this
    sample_list = read_sample(sample_file)
    if two_group_mode:
        if second_sample_file == "":
            print("Second sample file is missing - aborting 3")
            sys.exit(3)
        else:
            second_sample_list = read_sample(second_sample_file)
            sample_list = sample_list + second_sample_list
    df.index = sample_list  # change index to sample

    # get psi info

    info_df = pd.DataFrame()
    line_count = 0
    info_list = []
    with open(from_gtf, 'r') as f:
        for line in f:
            if line_count == 0:  # header
                line_count += 1
                info_index = line.rstrip().split('\t')[1:]
                continue
            elements = line.rstrip().split('\t')
            event_ID = str(elements[0])
            event_info = elements[1:]
            event_info_clean = [x.strip('\"') for x in event_info]  # genes and transcript ID has quotes - remove those
            new_info_df = pd.DataFrame({event_ID: event_info_clean})
            info_list.append(new_info_df)
            line_count += 1

    info_df = pd.concat(info_list, axis=1)

    info_df.index = info_index  # put name in
    # if event_type == "SE" or event_type == "MXE" or event_type == "RI":
    # 	info_df = info_df.drop(["upstreamES", "downstreamEE"]) # drop irrelevent info

    out_matrix = pd.concat([info_df, df], join='inner').T  # merge name and transpose

    print((info_df.index, df.index))

    # split by chromosome
    if split_by_chromosome:
        out_matrix.sort_values(by=['chr'], inplace=True)
        out_matrix.set_index(keys=['chr'], inplace=True, drop=False)
        chromosomes = out_matrix['chr'].unique().tolist()

        out_file_list = []
        for chromosome in chromosomes:
            # print chromosome
            chromosome_df = out_matrix.loc[out_matrix.chr == chromosome]
            num_chromosome_specific_event = chromosome_df.shape[0]
            print(("%i events are in %s" % (num_chromosome_specific_event, chromosome)))

            out_file_chromosome = "%s/merged_matrix.%s.%s_split.%s.txt" % (output_dir, event_type, tissue, chromosome)
            chromosome_df.to_csv(out_file_chromosome, na_rep='NaN', index=False, index_label=False)

            out_file_list.append(out_file_chromosome)

        out_file = ",".join(out_file_list)



    else:
        if chromosome == "":
            out_file = "%s/merged_matrix.%s.%s.txt" % (output_dir, event_type, tissue)
        else:
            out_file = "%s/merged_matrix.%s.%s_%s.txt" % (output_dir, event_type, tissue, chromosome)

        out_matrix.to_csv(out_file, na_rep='NaN', index=False, index_label=False)  # remove index for easier merge in R

    out_shape = out_matrix.shape
    print(("Matrix from %i samples and %i events " % (out_shape[1] - 8, out_shape[0])))
    if NA_filter:
        print(("%i events filtered out for missing more than %f percent of samples" % (num_below_NA, NA_cutoff * 100.0)))
    if coverage_filter:
        print(("%i events filtered out for low mean coverage below %f " % (num_below_coverage, coverage_cutoff)))

    # index test

    # new_index = out_matrix[:8].apply(",".join, axis = 0)

    # print new_index.shape
    return out_file


def test():
    from_gtf = "./test_sample/from_gtf.test.txt"
    jc_file = "./test_sample/jc.test.txt"
    sample_file = "./test_sample/samples.test.txt"
    tissue = "Blood"
    outdir = "./test_sample/"
    print((generate_psi_matrix(from_gtf, jc_file, sample_file, outdir, "SE", tissue, NA_filter=True,
                              split_by_chromosome=True)))


if __name__ == "__main__":

    # test()
    # sys.exit(4)

    script_dir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(description="PSI Matrix Computation Script for rMATS Output")

    input_args = parser.add_argument_group("Input Arguments")

    input_args.add_argument("input_dir", help="The dir where the rMATS results are stored.")
    input_args.add_argument("output_dir", help="The dir where the PSI matrix will be stored.")

    # optional arguments
    input_args.add_argument("--tissue", help="The tissue name - if not given, it will be selected automatically",
                            default="")
    input_args.add_argument("--event", help="The splicing event - default = Skipped Exon (SE)", default="SE")
    input_args.add_argument("--sample", help="Sample file if different", default="")
    input_args.add_argument("--two_group", help="If two groups are used - default = one group", action="store_true",
                            default=False)
    input_args.add_argument("--second_sample_file", help="Second sample file if you use two_group mode", default="")
    input_args.add_argument("--debug", help="Debug mode - do not use it", action="store_true")
    input_args.add_argument("--chromosome", help="run by chromosome, under development ", action="store", default="")

    # filter arguments
    filter_args = parser.add_argument_group("Filter Arguments")
    filter_args.add_argument("--na_filter", help="Apply NA Filter for the new matrix", action="store_true")
    filter_args.add_argument("--na_cutoff",
                             help="NA filter cutoff - any events missing more than the cutoff is removed - default = 0.9",
                             default=0.90)
    filter_args.add_argument("--coverage_filter", help="Apply mean cutoff filter for the new matrix",
                             action="store_true")
    filter_args.add_argument("--coverage_cutoff",
                             help="Coverage filter cutoff - anything with mean coverage below this threshold is removed - default = 10",
                             default=10.0)

    output_args = parser.add_argument_group("Output Arguments")
    output_args.add_argument("--split_chromosome", help="Split the output matrix by chromosome", action="store_true")
    output_args.add_argument("--event_annotation", help="Make event annotation file", action="store_true")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.exists(output_dir):  # make analysis dir if it does not exist
        os.makedirs(output_dir)

    if args.tissue == "":
        tissue = input_dir.strip('/').split('/')[-1].split('_rMATS_post')[0]
    else:
        tissue = args.tissue

    event = args.event

    if event not in ["SE", "RI", "MXE", "A3SS", "A5SS"]:
        print(("The input event - %s -  is not a valid event type - please check" % event))
        sys.exit(3)
    # TODO - add function for automatic run for all 5

    # select the files
    if args.chromosome == "":
        from_gtf_file = "%s/fromGTF.%s.txt" % (input_dir, event)
        JC_file = "%s/JC.raw.input.%s.txt" % (input_dir, event)
    else:  # if chromosome mode
        chromosome = args.chromosome
        from_gtf_file = "%s/fromGTF.%s_%s.txt" % (input_dir, event, chromosome)
        JC_file = "%s/JC.raw.input.%s_%s.txt" % (input_dir, event, chromosome)

    if args.sample == "":
        sample_file = "%s/sample.txt" % (input_dir)
    else:  # if sample file is given
        sample_file = args.sample

    if args.debug:
        psi_matrix = generate_psi_matrix_dev(from_gtf_file, JC_file, sample_file, output_dir, tissue, event,
                                             two_group_mode=args.two_group)
    else:
        psi_matrix = generate_psi_matrix(from_gtf_file, JC_file, sample_file,
                                         output_dir, tissue, event,
                                         two_group_mode=args.two_group, second_sample_file=args.second_sample_file,
                                         NA_filter=args.na_filter, NA_cutoff=args.na_cutoff,
                                         coverage_filter=args.coverage_filter, coverage_cutoff=args.coverage_cutoff,
                                         split_by_chromosome=args.split_chromosome, chromosome=args.chromosome)

    print(("PSI matrix is generated - %s " % psi_matrix))

    # generate annotation file for downstream analysis
    print("Generating annotation file")
    if args.event_annotation:
        anno_file = "%s/fromGTF.%s.txt" % (input_dir, event)
        noveljc_anno_file = "%s/fromGTF.novelJunction.%s.txt" % (input_dir, event)
        novelss_anno_file = "%s/fromGTF.novelSpliceSite.%s.txt" % (input_dir, event)

        print("Reading in the from GTF file")
        anno_df = pd.read_csv(anno_file, sep='\t', index_col=0)
        exon_list = [":".join([str(y) for y in x] + [event]) for index, x in anno_df.iterrows()]

        # read other ones
        noveljc_anno_df = pd.read_csv(noveljc_anno_file, sep='\t', index_col=0)
        novelss_anno_df = pd.read_csv(novelss_anno_file, sep='\t', index_col=0)

        # load index for annotation
        noveljc_id = set(noveljc_anno_df.index)
        novelss_id = set(novelss_anno_df.index)

        print("Annotating events")
        # generate annotation
        event_annotation = ["Novel_JC" if x in noveljc_id else "Novel_SS" if x in novelss_id else "Known_JC"
                            for x in anno_df.index]

        print("Joining the events")
        # make annotation df
        annotated_exon_df = pd.DataFrame([exon_list, event_annotation]).transpose()
        annotated_exon_df.columns = ['exon', 'annotation']
        annotation_out_file = "%s/AS_annotation.%s.%s.txt" % (output_dir, event, tissue)
        annotated_exon_df.to_csv(annotation_out_file, index=False, index_label=False)

        print("Output file stored in %s" % annotation_out_file)
