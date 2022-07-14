import pandas as pd
from glob import glob
import os
import argparse

#%%



#%%

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(description="")

    input_args = parser.add_argument_group("Input Arguments")
    input_args.add_argument("input_dir", help="Dir with kallisto run. Should contain dirs ending with *_quant")
    input_args.add_argument("output_dir", help="Where the merged and converted count would be stored")
    input_args.add_argument("study", help="study name")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    study = args.study

    # read kallisto conversion table
    kallisto_conversion_table_file = "/home/yangt3/xinglab/references/kallisto/homo_sapiens/transcripts_to_genes.txt"
    kallisto_conversion_df = pd.read_csv(kallisto_conversion_table_file, sep='\t', header=None,
                                         names=['transcript_ID','gene_ID', 'gene_symbol'])

    # input_kallisto_file = "/Users/harryyang/gtex_sqtl/expression_test/kallisto.TCGA-AB-3011-03A.tsv"
    sample_exp_df_list = []
    # input_dir = "/home/yangt3/xinglab/hotspot/expression/count/"

    input_kallisto_file_list = glob("%s/*/abundance.tsv" % input_dir)

    for input_kallisto_file in input_kallisto_file_list:
        # extract sample name
        sample = input_kallisto_file.strip().split('/')[-2].split('_quant')[0]
        print("Processing sample - %s" % sample)
        # read in kallisto file
        input_kallisto_df = pd.read_csv(input_kallisto_file, sep='\t')
        # add gene ID and symbol by merging
        annotated_kallisto_df = input_kallisto_df.set_index('target_id').join(kallisto_conversion_df.set_index('transcript_ID'))
        # merge by gene names
        gene_exp_df_list = []
        for gene_id, exp_block in annotated_kallisto_df.groupby('gene_ID'):
            exp_value = sum(exp_block['tpm'])
            gene_exp_df = pd.Series({'geneID': gene_id,
                                     'geneSymbol': exp_block['gene_symbol'][0],
                                     'TPM': exp_value})
            gene_exp_df_list.append(gene_exp_df)
        # merge the dataframes
        sample_exp_df = pd.concat(gene_exp_df_list, axis = 1).transpose()
        # set index - gene id and gene symbol is unique so it need to be set as index
        sample_exp_df.set_index(['geneID', 'geneSymbol'], inplace = True)
        # change the column as sample name so it can be utilzied for merging
        sample_exp_df.columns = [sample]
        sample_exp_df_list.append(sample_exp_df)

    #%%
    # merge the sample
    merged_df = pd.concat(sample_exp_df_list, axis = 1)

    # write the file
    output_file = "%s/%s.merged.exp.kallisto.txt" % (output_dir, study)
    merged_df.to_csv(output_file, sep='\t')
    print("Merged expression stored in - %s" % output_file)




