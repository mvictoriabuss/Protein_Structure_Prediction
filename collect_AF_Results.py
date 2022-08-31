#!/usr/bin/python

import argparse
import json
import os
import pandas as pd
import sys


def parse_arguments():
    """
    Short script to organize Alphafold results from a given set of proteins. Returns a summary score file used by ManageFoldseek.py
    * AlphaFold(AF) results root directory (AF_results_dir) and destination CSV (outfile) are required
    """
    parser = argparse.ArgumentParser(prog='Collect_AF_Results', description='Collects results from running of Alphafold to output\
                                      a summary CSV file.',
                                     usage='%(prog)s AF_results_dir, outfile')
    parser.add_argument('AF_results_dir', type=str, help='Directory containing the individual AF dirs of each protein')
    parser.add_argument('outfile', type=str, help='Path and name of csv output file')
        
    try:
        return parser.parse_args()
    except SystemExit:
        sys.exit(1)

def run():
    args = parse_arguments()
    
    #  Setup an Pandas dataframe to collect the data for each result. 
    df = pd.DataFrame(columns=['protein','best model', 'model score'])
    
    results_list = os.listdir(args.AF_results_dir)  # Results are in individual subdirectories, named by protein ID
    # Iterate thorugh results to access each result file, parsing the JSON to retrieve just the model scores 'plddts'
    for r in results_list:
        ranking_file = args.AF_results_dir + '/' + r + '/' + r + '.aa/ranking_debug.json'

        if os.path.isfile(ranking_file):
            with open(ranking_file, 'r') as file:
                model_scores = json.load(file)['plddts']
                
                #dict > list, and sort by the score which is the second element
                model_scores_list = sorted(model_scores.items(), key=lambda x:x[1])
                
                #Best model is the last in the list, we want the model # and the corresponding score
                top_model, top_score = model_scores_list[-1][0], model_scores_list[-1][1]
                result = {'protein': r, 'best model': top_model, 'model score': top_score}
                df.loc[len(df.index)] = result
        else:
            result = {'protein': r, 'best model': 'missing', 'model score': 'missing'}
            df.loc[len(df.index)] = result

    df = df.set_index('protein')
    df.to_csv(args.outfile)

if __name__ == '__main__':
    run()
    