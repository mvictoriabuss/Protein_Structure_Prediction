#!/usr/bin/python

import csv
import os
import pandas as pd
import shutil
import subprocess
import sys


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

def main(AF_results_dir, score_file, best_pdb_dir, foldseek_db, FS_output_dir, AFcutoff=70, FScutoff='1.000E-5', threads=24):
	make_saved_pdb_dir(best_pdb_dir)
	make_output_dir(FS_output_dir)
	df = read_scores(score_file, int(AFcutoff))
	copy_pdb(df, AF_results_dir, best_pdb_dir)
	seeklist = list_PDBs(best_pdb_dir)

	options = f' --threads {threads} --format-output target,evalue,tlen,qlen,qstart,qend,tcov,qcov,taxname -e {FScutoff}'
	run_Foldseek(seeklist, best_pdb_dir, foldseek_db, FS_output_dir, options)

	dpr = make_fs_summary_CSV(FS_output_dir, 'FS_summary.csv')  # Dataframe PDB Results
	rerun_list = make_rerun_list(dpr)
	
	foldseek_db = foldseek_db.replace('pdb', 'afdb')
	FS_output_dir2 = FS_output_dir + '2nd_DB/'
	make_output_dir(FS_output_dir2)
	run_Foldseek(seeklist, best_pdb_dir, foldseek_db, FS_output_dir2, options)
	_ = make_fs_summary_CSV(FS_output_dir2, 'FS_summary_afdb.csv')

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print('Usage: AF_results_dir, score_file, best_pdb_dir, foldseek_db, output_dir')
		raise SyntaxError("Insufficient arguments.")
	else:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
