import pandas as pd
import pybedtools
import gzip
from glob import glob
import os
import argparse
import sys

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(description="")

    input_args = parser.add_argument_group("Input Arguments")
    input_args.add_argument("input_psi_file",
                            help="Filtered PSI matrix from hotspot analysis step")
    input_args.add_argument('output_dir',
                            help="directory to save maf files")
    input_args.add_argument('cancer_type',
                            help="cancer type from TCGA - e.g. TCGA-BRCA")
    input_args.add_argument('event_type',
                            help="alternative splicing event type - eg. SE")
    input_args.add_argument('--encode_bed',
                            help="Dir with Peak BED files from ENCODE",
                            default='/home/yangt3/xinglab/ENCODE/bed/')
    input_args.add_argument('--underscore',
                            help="Change process for Exon formatted with underscore",
                            action="store_true")

    args = parser.parse_args()

    # # local debug
    # encode_bed_dir = "/Users/harryyang/research/Xing Lab/TumorAnalyisFramework/encode/peaks/"
    # output_dir = "/Users/harryyang/research/Xing Lab/TumorAnalyisFramework/encode/test_3/"
    # # input_psi_file = '/Users/harryyang/research/Xing Lab/TumorAnalyisFramework/encode/test.csv'
    # input_psi_file = '/Users/harryyang/research/Xing Lab/TumorAnalyisFramework/data/YangPan_PNAS/myc_cellline_deltapsi.csv'
    # cancer_type = "TCGA-SKCM"
    # event_type = "SE"
    # underscore_mode = True

    input_psi_file = args.input_psi_file
    output_dir = args.output_dir
    cancer_type = args.cancer_type
    event_type = args.event_type
    encode_bed_dir = args.encode_bed
    underscore_mode = args.underscore

    output_intersect_dir = "%s/intersect/" % output_dir
    # initialization
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(output_intersect_dir):
        os.mkdir(output_intersect_dir)

    # Process exon list
    intron_padding = 300
    exon_padding = 50

    # import exons into bed
    print("Reading in PSI table for AS events")
    psi_df = pd.read_csv(input_psi_file)
    if not underscore_mode:
        exon_list = [":".join([str(y) for y in x[:10]] + [event_type]) for index, x in psi_df.iterrows()]
    else:
        exon_list = psi_df.iloc[:, 0].tolist()

    # generate exonic regions
    print("Generating regions of interest for AS events")
    gff_list = []
    event_count = 0
    interval_bed_file = "%s/exon_region.%s.gff.gz" % (output_dir, cancer_type)
    source = "TumorAnalysisFramework"

    with gzip.open(interval_bed_file, 'wb') as f:
        for exon in exon_list:
            if event_type == "SE":
                if not underscore_mode:
                    exon_elements = exon.split(":")

                    chromosome = exon_elements[2]
                    strand = exon_elements[3]

                    # define ss junctions and add entries
                    if strand == "+":
                        upstream_5ss = int(exon_elements[7])
                        cassette_3ss = int(exon_elements[4])
                        cassette_5ss = int(exon_elements[5])
                        dnstream_3ss = int(exon_elements[8])
                    elif strand == "-":
                        upstream_5ss = int(exon_elements[8])
                        cassette_3ss = int(exon_elements[5])
                        cassette_5ss = int(exon_elements[4])
                        dnstream_3ss = int(exon_elements[7])
                elif underscore_mode:
                    exon_elements = exon.split("_")

                    chromosome = exon_elements[1]
                    strand = exon_elements[2]

                    # define ss junctions and add entries
                    if strand == "+":
                        upstream_5ss = int(exon_elements[5])
                        cassette_3ss = int(exon_elements[3])
                        cassette_5ss = int(exon_elements[4])
                        dnstream_3ss = int(exon_elements[6])
                    elif strand == "-":
                        upstream_5ss = int(exon_elements[6])
                        cassette_3ss = int(exon_elements[4])
                        cassette_5ss = int(exon_elements[3])
                        dnstream_3ss = int(exon_elements[5])

                if strand == "+":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - exon_padding, upstream_5ss, ".", strand, 0,
                         'feature=Upstream_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + intron_padding, ".", strand, 0,
                         'feature=Upstream_5ss_intron'])
                    cassette_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_3ss - intron_padding, cassette_3ss, ".", strand, 0,
                         'feature=Cassette_3ss_intron'])
                    cassette_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_3ss, cassette_3ss + exon_padding, ".", strand, 0,
                         'feature=Cassette_3ss_exon'])
                    cassette_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_5ss - exon_padding, cassette_5ss, ".", strand, 0,
                         'feature=Cassette_5ss_exon'])
                    cassette_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_5ss, cassette_5ss + intron_padding, ".", strand, 0,
                         'feature=Cassette_5ss_intron'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - intron_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Dnstream_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + exon_padding, ".", strand, 0,
                         'feature=Dnstream_3ss_exon'])

                elif strand == "-":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + exon_padding, ".", strand, 0,
                         'feature=Upstream_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - intron_padding, upstream_5ss, ".", strand, 0,
                         'feature=Upstream_5ss_intron'])
                    cassette_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_3ss, cassette_3ss + intron_padding, ".", strand, 0,
                         'feature=Cassette_3ss_intron'])
                    cassette_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_3ss - exon_padding, cassette_3ss, ".", strand, 0,
                         'feature=Cassette_3ss_exon'])
                    cassette_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_5ss, cassette_5ss + exon_padding, ".", strand, 0,
                         'feature=Cassette_5ss_exon'])
                    cassette_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, cassette_5ss - intron_padding, cassette_5ss, ".", strand, 0,
                         'feature=Cassette_5ss_intron'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + intron_padding, ".", strand, 0,
                         'feature=Dnstream_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - exon_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Dnstream_3ss_exon'])

                # add relevant entry as string to import to make a new bed
                combined_gff_entry = "%s%s%s%s%s%s%s%s" % (upstream_5ss_exon_interval,
                                                           upstream_5ss_intr_interval,
                                                           cassette_3ss_intr_interval,
                                                           cassette_3ss_exon_interval,
                                                           cassette_5ss_exon_interval,
                                                           cassette_5ss_intr_interval,
                                                           dnstream_3ss_intr_interval,
                                                           dnstream_3ss_exon_interval)

            elif event_type == "A3SS":
                print(exon)
                if not underscore_mode:
                    exon_elements = exon.split(":")

                    chromosome = exon_elements[2]
                    strand = exon_elements[3]

                    # define ss junctions and add entries
                    if strand == "+":
                        upstream_5ss = int(exon_elements[9])
                        alternat_3ss = int(exon_elements[4])
                        dnstream_3ss = int(exon_elements[6])
                    elif strand == "-":
                        upstream_5ss = int(exon_elements[8])
                        alternat_3ss = int(exon_elements[5])
                        dnstream_3ss = int(exon_elements[7])
                elif underscore_mode:
                    print("Underscore mod not supported for A3SS")
                    sys.exit(1)

                if strand == "+":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - exon_padding, upstream_5ss, ".", strand, 0,
                         'feature=Upstream_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + intron_padding, ".", strand, 0,
                         'feature=Upstream_5ss_intron'])
                    alternat_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_3ss - intron_padding, alternat_3ss, ".", strand, 0,
                         'feature=Long_3ss_intron'])
                    alternat_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_3ss, alternat_3ss + exon_padding, ".", strand, 0,
                         'feature=Long_3ss_exon'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - intron_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Short_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + exon_padding, ".", strand, 0,
                         'feature=Short_3ss_exon'])

                elif strand == "-":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + exon_padding, ".", strand, 0,
                         'feature=Upstream_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - intron_padding, upstream_5ss, ".", strand, 0,
                         'feature=Upstream_5ss_intron'])
                    alternat_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_3ss, alternat_3ss + intron_padding, ".", strand, 0,
                         'feature=Long_3ss_intron'])
                    alternat_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_3ss - exon_padding, alternat_3ss, ".", strand, 0,
                         'feature=Long_3ss_exon'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + intron_padding, ".", strand, 0,
                         'feature=Short_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - exon_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Short_3ss_exon'])

                # add relevant entry as string to import to make a new bed
                combined_gff_entry = "%s%s%s%s%s%s" % (upstream_5ss_exon_interval,
                                                       upstream_5ss_intr_interval,
                                                       alternat_3ss_intr_interval,
                                                       alternat_3ss_exon_interval,
                                                       dnstream_3ss_intr_interval,
                                                       dnstream_3ss_exon_interval)

            elif event_type == "A5SS":
                print(exon)
                if not underscore_mode:
                    exon_elements = exon.split(":")

                    chromosome = exon_elements[2]
                    strand = exon_elements[3]

                    # define ss junctions and add entries
                    if strand == "+":
                        upstream_5ss = int(exon_elements[7])
                        alternat_5ss = int(exon_elements[5])
                        dnstream_3ss = int(exon_elements[8])
                    elif strand == "-":
                        upstream_5ss = int(exon_elements[6])
                        alternat_5ss = int(exon_elements[4])
                        dnstream_3ss = int(exon_elements[9])
                elif underscore_mode:
                    print("Underscore mod not supported for A5SS")
                    sys.exit(1)

                if strand == "+":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - exon_padding, upstream_5ss, ".", strand, 0,
                         'feature=Short_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + intron_padding, ".", strand, 0,
                         'feature=Short_5ss_intron'])
                    alternat_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_5ss - exon_padding, alternat_5ss, ".", strand, 0,
                         'feature=Long_5ss_exon'])
                    alternat_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_5ss, alternat_5ss + intron_padding, ".", strand, 0,
                         'feature=Long_5ss_intron'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - intron_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Dnstream_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + exon_padding, ".", strand, 0,
                         'feature=Dnstream_3ss_exon'])

                elif strand == "-":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + exon_padding, ".", strand, 0,
                         'feature=Short_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - intron_padding, upstream_5ss, ".", strand, 0,
                         'feature=Short_5ss_intron'])
                    alternat_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_5ss, alternat_5ss + exon_padding, ".", strand, 0,
                         'feature=Long_5ss_exon'])
                    alternat_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, alternat_5ss - intron_padding, alternat_5ss, ".", strand, 0,
                         'feature=Long_5ss_intron'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + intron_padding, ".", strand, 0,
                         'feature=Dnstream_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - exon_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Dnstream_3ss_exon'])

                # add relevant entry as string to import to make a new bed
                combined_gff_entry = "%s%s%s%s%s%s" % (upstream_5ss_exon_interval,
                                                       upstream_5ss_intr_interval,
                                                       alternat_5ss_exon_interval,
                                                       alternat_5ss_intr_interval,
                                                       dnstream_3ss_intr_interval,
                                                       dnstream_3ss_exon_interval)
            elif event_type == "RI":
                print(exon)
                if not underscore_mode:
                    exon_elements = exon.split(":")

                    chromosome = exon_elements[2]
                    strand = exon_elements[3]

                    # define ss junctions and add entries
                    if strand == "+":
                        upstream_5ss = int(exon_elements[4])
                        dnstream_3ss = int(exon_elements[5])
                    elif strand == "-":
                        upstream_5ss = int(exon_elements[5])
                        dnstream_3ss = int(exon_elements[4])
                elif underscore_mode:
                    print("Underscore mod not supported for RI")
                    sys.exit(1)

                if strand == "+":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - exon_padding, upstream_5ss, ".", strand, 0,
                         'feature=Upstream_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + intron_padding, ".", strand, 0,
                         'feature=Upstream_5ss_intron'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - intron_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Dnstream_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + exon_padding, ".", strand, 0,
                         'feature=Dnstream_3ss_exon'])

                elif strand == "-":
                    # make new interval files
                    upstream_5ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss, upstream_5ss + exon_padding, ".", strand, 0,
                         'feature=Upstream_5ss_exon'])
                    upstream_5ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, upstream_5ss - intron_padding, upstream_5ss, ".", strand, 0,
                         'feature=Upstream_5ss_intron'])
                    dnstream_3ss_intr_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss, dnstream_3ss + intron_padding, ".", strand, 0,
                         'feature=Dnstream_3ss_intron'])
                    dnstream_3ss_exon_interval = pybedtools.bedtool.create_interval_from_list(
                        [chromosome, source, exon, dnstream_3ss - exon_padding, dnstream_3ss, ".", strand, 0,
                         'feature=Dnstream_3ss_exon'])

                # add relevant entry as string to import to make a new bed
                combined_gff_entry = "%s%s%s%s" % (upstream_5ss_exon_interval,
                                                   upstream_5ss_intr_interval,
                                                   dnstream_3ss_intr_interval,
                                                   dnstream_3ss_exon_interval)
            exon_gff = pybedtools.bedtool.BedTool(combined_gff_entry, from_string=True)

            # write them
            f.write(str(exon_gff).encode())
            event_count += 1
            if event_count % 1000 == 0:
                print("Procesesed %i events for interval generation" % event_count)

    # intersect for each RBP peaks
    print("Intersecting AS regions with CLIP peaks")
    intersect_count_list = []
    input_clip_peak_file_list = glob("%s/*.bed.gz" % encode_bed_dir)
    for input_clip_peak_file in input_clip_peak_file_list:
        input_clip_peak_bed = pybedtools.BedTool(input_clip_peak_file)

        # infer CLIP data
        study = input_clip_peak_file.split('/')[-1].split('.')[0]
        rbp = input_clip_peak_file.split('/')[-1].split('.')[1]

        print("Processing %s - %s" % (study, rbp))

        # read in the merged file
        # writing and reading is easier as bedtools merge can only run one entry at a time
        # - simplification via reduce takes too much time
        merged_gff = pybedtools.bedtool.BedTool(interval_bed_file)
        gff_clip_intersect = merged_gff.intersect(input_clip_peak_bed)
        intersect_output_file = "%s/%s.%s.%s.gff.gz" % (
            output_intersect_dir, cancer_type, study, rbp)  # TODO - fix this

        with gzip.open(intersect_output_file, 'wb') as f:
            f.write(str(gff_clip_intersect).encode())

        # import it back to pandas dataframe - these IO steps are easier to port between data formats
        try:
            intersect_df = pd.read_csv(intersect_output_file, sep='\t', header=None)
        except:  # empty ones can happen if there is no overlap
            continue
        intersect_df.columns = ['chr', 'study', 'exon', 'start', 'end', 'score', 'strand', 'orf', 'feature']

        # trim feature
        intersect_df['region'] = [x.split('feature=')[1] for x in intersect_df['feature']]

        # drop duplicates just in case
        intersect_df = intersect_df.drop_duplicates()
        intersect_df.names = 'count'

        # generate count table
        intersect_count = intersect_df.value_counts(['exon', 'region']).reset_index(drop=False)

        # this analysis may be deferred to downstream
        # count_table = intersect_count.pivot(index = 'exon', columns = 'region', values = 'count')

        intersect_count['study'] = study
        intersect_count['RBP'] = rbp
        intersect_count['event_type'] = event_type

        intersect_count_list.append(intersect_count)

    # merge the counts
    print("Merging intersection results")
    intersect_count_df = pd.concat(intersect_count_list).reset_index(drop=True)
    intersect_count_df.columns = ['exon', 'region', 'count', 'study', 'RBP', 'event_type']
    output_file = '%s/%s.CLIP_overlap.%s.tsv' % (output_dir, cancer_type, event_type)
    intersect_count_df.to_csv(output_file, sep='\t')

    print("The output file is stored in %s" % output_file)
