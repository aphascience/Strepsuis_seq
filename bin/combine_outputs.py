#! /usr/bin/env python3

# Combine outputs from all four srst2 runs to generate a single results file in csv format

import pandas as pd
import argparse
import getpass
from datetime import date


def combineData(recNTable, MLSTTable, serotypeTable, virulenceTable, verifyCSV):

    # Get info for logging
    date_out = date.today().strftime('%d%b%y')
    user = getpass.getuser()

    # read recN data
    recN_df = pd.read_table(recNTable, sep='\t')
    recN_df['Sample'] = recN_df['Sample'].astype(object)
    recN_df['recN-Pos'].fillna('Not Ssuis', inplace=True)

    # read MLST data
    MLST_df = pd.read_table(MLSTTable, sep='\t')
    MLST_df[['Sample', 'ST', 'aroA', 'cpn60', 'dpr', 'gki', 'mutS', 'recA', 'thrA']] = \
        MLST_df[['Sample', 'ST', 'aroA', 'cpn60', 'dpr', 'gki', 'mutS', 'recA', 'thrA']].astype(object)

    # read serotype data
    serotype_df = pd.read_table(serotypeTable, sep='\t')
    serotype_df['Sample'] = serotype_df['Sample'].astype(object)
    serotype_df.rename({'ST': 'Serotype'}, axis=1, inplace=True)
    serotype_df['Serotype'] = serotype_df['Serotype'].str.replace(r'\D', '', regex=True)

    # read serotype verification
    verify_df = pd.read_csv(verifyCSV)
    verify_df[['Sample', 'Pos483', 'RefF', 'RefR', 'AltF', 'AltR']] = \
        verify_df[['Sample', 'Pos483', 'RefF', 'RefR', 'AltF', 'AltR']].astype(object)

    # update serotypes based on verification
    verified_serotype_df = pd.merge(serotype_df, verify_df, on='Sample', how='left')
    # Quote from https://rdcu.be/d56Ee:
    # ## "The analysis revealed that all serotype 2 and all serotype 14 strains had a G nucleotide at position 483 of
    # ## the cpsK gene, while all serotype 1 and all serotype 1/2 strains contained either a C or T at that nucleotide
    # ## position."
    # if prelim serotype is 1 and base 483 is G, then Serotype is 14
    # if prelim serotype is 1 and base 483 is C or T, then Serotype is 1
    # if prelim serotype is 2 and base 483 is G, then Serotype is 2
    # if prelim serotype is 2 and base 483 is C or T, then Serotype is 1/2

    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '1') & (verified_serotype_df['Pos483'] == 'G'),
                             'Serotype'] = '14'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '1') & (verified_serotype_df['Pos483'] == 'C'),
                             'Serotype'] = '1'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '1') & (verified_serotype_df['Pos483'] == 'T'),
                             'Serotype'] = '1'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '2') & (verified_serotype_df['Pos483'] == 'G'),
                             'Serotype'] = '2'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '2') & (verified_serotype_df['Pos483'] == 'C'),
                             'Serotype'] = '1/2'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '2') & (verified_serotype_df['Pos483'] == 'T'),
                             'Serotype'] = '1/2'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] != '1') & (verified_serotype_df['Serotype'] != '2'),
                             'Serotype'] = verified_serotype_df['Serotype']

    # read virulence data
    virulence_df = pd.read_table(virulenceTable, sep='\t', names=list(range(4)), skiprows=1)
    virulence_df.rename({0: 'Sample'}, axis=1, inplace=True)
    virulence_df.fillna('-', inplace=True)
    virulence_df['epf'] = virulence_df[1].str.contains('epf')
    virulence_df['mrp'] = (virulence_df[1].str.contains('mrp')) | (virulence_df[2].str.contains('mrp')) | \
        (virulence_df[3].str.contains('mrp'))
    virulence_df['sly'] = (virulence_df[1].str.contains('sly')) | (virulence_df[2].str.contains('sly')) | \
        (virulence_df[3].str.contains('sly'))
    virulence_df = virulence_df[['Sample', 'epf', 'mrp', 'sly']]

    # Merge dataframes
    sero_mlst_df = pd.merge(verified_serotype_df, MLST_df, on='Sample', how='left')

    finalout_df = pd.merge(pd.merge(recN_df, sero_mlst_df, on='Sample', how='left'),
                           virulence_df, on='Sample', how='outer')
    finalout_df.set_index('Sample', inplace=True)
    finalout_df.sort_index(inplace=True)

    # Write to csv
    finalout_df.to_csv("FinalOut_{}.csv".format(date_out))

    # Append log info
    with open("FinalOut_{}.csv".format(date_out), "a") as outFile:
        outFile.write("# Operator: " + user)
        outFile.close


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('recNTable', help='............')
    parser.add_argument('MLSTTable', help='............')
    parser.add_argument('serotypeTable', help='...............')
    parser.add_argument('virulenceTable', help='................')
    parser.add_argument('verifyCSV', help='................')
    # parser.add_argument('commitId', help='Nextflow capture of git commit')

    args = parser.parse_args()

    combineData(**vars(args))
