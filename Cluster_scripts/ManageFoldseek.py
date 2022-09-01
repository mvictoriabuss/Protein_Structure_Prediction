#!/usr/bin/python

import argparse
import csv
import os
import pandas as pd
import shutil
import subprocess
import sys

def parse_args():
	cwd = os.getcwd()
	parser = argparse.ArgumentParser(prog='Manage FoldSeek', description='Runs FoldSeek on predicted AF structures using structural DBs ',
                                     usage='%(prog)s AF_results_dir, FoldSeek_db, Score_File')
	parser.add_argument('AF_Results_Dir', type=str, help='Directory contains AF results')
	parser.add_argument('FoldSeek_DB', type=str, help='FoldSeek structural Database')
	parser.add_argument('Score_File', type=str, help="CSV from 'Collect_AF_Results.py', listing best models and scores")
	parser.add_argument('-o', '--outdir', help="Output files directory. Default=output_FoldSeek/", type=str,
                        default=f'{cwd}/output_FoldSeek')
	parser.add_argument('-p', '--top_models', help='Location to which best AF models, PDB format, will be copied: Default=top_AF_models/', 
                        type=str, default=f'{cwd}/top_AF_models')
	parser.add_argument('-t','--AF_score_threshold', type=int, default=80,
						help='Only AF models scores equal or greater than the threshold will be considered. Default: 80')
	parser.add_argument('-c', '--FS_score_cutoff', type=str, default='1.000E-5',
						help="As a string, the cutoff for return results from FoldSeek. Default: '1.000E-5'")
	parser.add_argument('-p', '--processors', type=int, default=12,
						help="Number of processors used for computation. Default: 12")
	
	try:
		args = parser.parse_args()
	except SystemExit:
		sys.exit(1)
		
	return args

def read_result(resultfile):
    with open(resultfile) as file:
        file_length = len(file.readlines())
    return file_length

def make_saved_pdb_dir(best_pdb_dir):
	if not os.path.exists(best_pdb_dir):
		try:
			os.mkdir(best_pdb_dir)
		except Exception as e:
			raise e
	else:
	    print(f'Best PDB directory, {best_pdb_dir}, already exists.\n Just sayin...')

def make_output_dir(output_dir):
	if not os.path.exists(output_dir):
		try:
			os.mkdir(output_dir)
		except Exception as e:
			raise e
	else:
	    print(f'Output directory, {output_dir} exists!\n Overwriting...')

def copy_pdb(df, AF_results_dir, best_pdb_dir):
	#  Copy the best model PDB files in a single directory, adding the protein name.
	for i in df.index:
		bestModel = df.loc[i]['best model']  # One of the five generated has the best score.
		score = int(round(df.loc[i]['model score']))
		pdb_file = AF_results_dir + '/' + i + '/' + i + '.aa/relaxed_' + bestModel + '.pdb'
		file_copy_rename = best_pdb_dir + '/' + i + '_modeled_' + str(score) + '.pdb'
		shutil.copy2(pdb_file, file_copy_rename)

def read_scores(score_file, AFcutoff):
	dx = pd.read_csv(score_file)
	dx = dx.set_index('protein')
	dx = dx[(dx['model score'] != 'missing')].copy()  # Remove any entries marked as missing. Make copy due to slicing issues.
	dx['model score'] = pd.to_numeric(dx['model score'])  # Make sure this is a float, not a str
	df = dx[(dx['model score'] >= AFcutoff)]
	return df

def list_PDBs(best_pdb_dir):
	dirlist = os.listdir(best_pdb_dir)
	seeklist = [i for i in dirlist if i.endswith('pdb')]
	return seeklist

def run_Foldseek(seeklist, best_pdb_dir, foldseek_db, output_dir, options):
	counter = 0
	for file in seeklist:
		protein_name = file.split('_modeled')[0]
		resultfile = output_dir + '/' + protein_name + '.m8'
		command = 'foldseek easy-search ' + best_pdb_dir + file + ' ' + foldseek_db + ' ' + resultfile + ' tmp' + options + ' >log.txt'
		process = subprocess.call(command, shell = True)
		file_length = read_result(resultfile)
		counter += 1
		print(f'Done with {protein_name}: {file_length} hits. {counter}/{len(seeklist)} now completed.')

def filter_results(results):
	"""Filter each FoldSeek result file from (make_fs_summary_CSV). Number of hits, empty results, and removal of
	hit to human are handled.

	Args:
		results (list[str]): The individuals results from a given PDB structure query result from Foldseek

	Returns:
		results[counter] (str): Selected result after filtering 
		number_of_hits (int): Number of hits in result file.
		note (str): Relevant information about reported hits.

	"""
	
	number_of_hits = len(results)
	counter = 0
    
    #If the file it's empty
	if number_of_hits == 0:
		number_of_hits = 'No hits on PDB'
		note = ''
		return False, number_of_hits, note

	#If it's not empty 
	else:
		#While loop to organise the results and look for Homo sapiens hits 
		while counter < len(results) and 'sapiens' in results[counter]:
			counter+=1

		#If file contains only Homo sapiens hits in PDB  
		if counter == len(results):
			note = "Only Homo Sapiens hits on PDB"
			return False, number_of_hits, note 

		#If the first hit isn't Homo sapiens
		elif counter == 0:
			note = '' 
			return results[counter], number_of_hits, note   

		#If n hits were Homo sapiens but then another specie is found
		else:
			note = 'Result n° ' + str(counter + 1) 
			return results[counter], number_of_hits, note 

def make_fs_summary_CSV(FS_output_dir, output_csv_file):
	'''Generate a CSV to summarize results. A seperate function will exclude some Eukaryotic species.
		Return a corresponding Pandas Dataframe.'''
	# open the cvs file
	with open(os.path.join(output_csv_file),'w', newline='', encoding='utf-8') as csvfile: 

		# Setup Output CSV, label columns
		csvwriter = csv.writer(csvfile) 
		fields = ['Sequence', 'Hits', '1st protein hit', '1st E-value', '1st Specie', 'Notes']
		csvwriter.writerow(fields) 

        #open the folder with the Foldseek results
		for sequence in os.listdir(FS_output_dir):
			with open(os.path.join(FS_output_dir, sequence)) as result_files:
				rows = result_files.readlines()        

                #Filter the results
				result, number_of_hits, note = filter_results(rows)  

                #If there are good hits on PDB --> Organize the results
				if result != False:
					column = result.rstrip('\n').split('\t')
					protein = column[0]
					e_value = column[1]
					specie = column[8]

                #if it's False --> no hits OR only Homo Sapiens hits on PDB
				else:
					protein = ""
					e_value = ""
					specie = ""
				
				rows = [[sequence.replace('.m8',''), number_of_hits, protein, e_value, specie, note]]

                # writing the data rows
				csvwriter.writerows(rows)
	df = pd.read_csv(output_csv_file)
	return df

def make_rerun_list(dpr):
	df = dpr[(~dpr['1st protein hit'].notna())]  # use new df so we can update the orignal later.
	rereun_list = list(df['Sequence'])
	print(rereun_list)
	return rereun_list

def run():
	args = parse_args()

	df = read_scores(args.Score_File, int(args.AF_score_threshold))
	
	make_saved_pdb_dir(args.top_models)
	copy_pdb(df, args.AF_Results_Dir, args.top_models)
	seeklist = list_PDBs(args.top_models)
	
	make_output_dir(args.outdir)
	foldseek_db = args.FoldSeek_DB
	options = f' --threads {args.processors} --format-output target,evalue,tlen,qlen,qstart,qend,tcov,qcov,taxname -e {args.FS_score_cutoff}'
	run_Foldseek(seeklist, args.top_models, foldseek_db, args.outdir, options)

	dpr = make_fs_summary_CSV(args.outdir, 'FS_summary_pdb.csv')  # Dataframe PDB Results
	rerun_list = make_rerun_list(dpr)
	
	foldseek_db = foldseek_db.replace('pdb', 'afdb')
	FS_output_dir2 = args.outdir + '2nd_DB/'
	make_output_dir(FS_output_dir2)
	run_Foldseek(seeklist, args.top_models, foldseek_db, FS_output_dir2, options)
	_ = make_fs_summary_CSV(FS_output_dir2, 'FS_summary_afdb.csv')

if __name__ == "__main__":
	run()
