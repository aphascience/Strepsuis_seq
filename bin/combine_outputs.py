#! /usr/bin/env python3

# Combine outputs from all four srst2 runs to generate a single results file in csv format

import pandas as pd
import argparse
from datetime import datetime, date


def combineData(recNTable, MLSTTable, serotypeTable, virulenceTable):
    # read recN data
    recN_df = pd.read_table(recNTable, sep='\t')
    recN_df['Sample'] = recN_df['Sample'].astype(object)
    recN_df['recN-Pos'].fillna('Not Ssuis', inplace=True)
    print(recN_df)

    # read MLST data
    MLST_df = pd.read_table(MLSTTable, sep='\t')
    MLST_df['Sample'] = MLST_df['Sample'].astype(object)

    # read serotype data
    serotype_df = pd.read_table(serotypeTable, sep='\t')
    serotype_df['Sample'] = serotype_df['Sample'].astype(object)
    print(serotype_df)

    # read virulence data
    virulence_df = pd.read_table(virulenceTable, sep='\t', names=list(range(4)), skiprows=1)
    virulence_df.rename({0 : 'Sample'}, axis=1, inplace=True)
    virulence_df.fillna('-', inplace=True)
    virulence_df['epf'] = virulence_df[1].str.contains('epf')
    virulence_df['mrp'] = (virulence_df[1].str.contains('mrp')) | (virulence_df[2].str.contains('mrp')) | (virulence_df[3].str.contains('mrp'))
    virulence_df['sly'] = (virulence_df[1].str.contains('sly')) | (virulence_df[2].str.contains('sly')) | (virulence_df[3].str.contains('sly'))
    virulence_df = virulence_df[['Sample', 'epf', 'mrp', 'sly']]
    print(virulence_df)

    # Merge dataframes
    sero_mlst_df = pd.merge(MLST_df, serotype_df, on='Sample', how='left')

    finalout_df = pd.merge(pd.merge(recN_df, virulence_df, on='Sample', how='left'), sero_mlst_df, on='Sample', how='outer')
    finalout_df.set_index('Sample', inplace=True)
    finalout_df.sort_index(inplace=True)
    print(finalout_df)

    #Write to csv
    finalout_df.to_csv("FinalOut.csv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('recNTable', help='............')
    parser.add_argument('MLSTTable', help='............')
    parser.add_argument('serotypeTable', help='...............')
    parser.add_argument('virulenceTable', help='................')
    #parser.add_argument('commitId', help='Nextflow capture of git commit')
    #parser.add_argument('--read_threshold', type=int, default=500, help='threshold for number of M.bovis reads')
    #parser.add_argument('--abundance_threshold', type=int, default=1, help='threshold for M.bovis abundance')

    args = parser.parse_args()

    combineData(**vars(args))
