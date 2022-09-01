#!/usr/bin/python

import argparse
import gffpandas.gffpandas as gp
import os
import pyfaidx
import re
import sys
from datetime import datetime


def parse_arguments():
    """
    Parser for commandline arguments:
    * FASTA sequence (a multi-FASTA seqFile) is required a the first positional argument.
    * If no targetsFile is specified, the '-a' must also be used and all ORFs will be extracted.
    * Setting an output directory for Bash & translated FASTAs, desides the default, is optional.
    """
    parser = argparse.ArgumentParser(prog='Extract_Translate_Batch', description='Select ORFs for translation and Batch files.',
                                     usage='%(prog)s seqFile, targetsFile, outdir (optional)')
    parser.add_argument('seqFile', type=str, help='Multi FASTA file with all ORFs from the Genome.')
    parser.add_argument('-t', '--targetsFile', type=str, help='Single Column CSV contain ORFs for extraction',
                        default=None)
    parser.add_argument('-o', '--outdir', help="Output files directory. Default=output_aa_batch/", type=str,
                        default='output_aa_batch')
    parser.add_argument('-l', '--AF_out', help='Set the desired destination dir that Batch files will use for AF structures', 
                        type=str, default=None)
    parser.add_argument('-f', '--FASTA_Format', help='Sets the (header) format expected from the FASTA seqFile. Default=none',
                        type=str, choices=['refseq', 'genbank', 'gff', 'none'], default='none')
    parser.add_argument('-s', '--skip_prompt', help='Suppress output directory exists prompt.', action='store_true')
    parser.add_argument('-a', '--all_ORFs', help='Extracts all ORFs in seqFile, ignoring targetsFile which is not required.',
                        action='store_true', default=False)
     
    try:
        args = parser.parse_args()
    except SystemExit:
        sys.exit(1)

    if not args.targetsFile and not args.all_ORFs:
        raise Exception('If no TargetsFile is specified, then the all_ORFs option is required.')
    
    return args

def handle_outdir(outdir_path, suppress_outdir_prompt) -> str:
    """
    Check if the output directory exists.
    Unless suppressed with -os argument, user will need to confirm overwriting to continue
    :return: path to output directory
    """
    yes_answers = ['Y', 'YES', 'MAKE IT SO', 'OUI']
    if os.path.isdir(outdir_path):
        if suppress_outdir_prompt:
            print('Output directory %s/ already exists, but away we go anyway!' % outdir_path)
        else:
            answer = input('Output directory %s/ already exists! Continue anyway [y/N]?' % outdir_path).upper()
            if answer in yes_answers:
                print('Okay. Files in %s/ will be overwritten.' % outdir_path)
            else:
                raise SystemExit("Not overwriting...aborting.")

    # Create directory and report success.
    if not os.path.isdir(outdir_path):
        try:
            os.mkdir(outdir_path)
        except OSError:
            print('There was a problem creating %s/' % outdir_path)
        else:
            print('Output directory %s/ successfully created.' % outdir_path)
    return outdir_path

def make_targets(target_file, ORFs) -> list:
    """Here we create the list of target ORFs that will be extrated.
    If the TargetsFile argument is not specified, it defaults to None. In this case all ORFs are returned as targets.

    Args:
        target_file (str): The list of ORFs to extract, or 'None' for all.
        ORFs (pyfaidx.Fasta): Searchable format generated from multi-FASTA seqeunce file.

    Returns:
        list: The list of ORFs that will be extracted.
    """
    if target_file:
        with open(target_file) as file:
            targets = file.readlines()

        targets = ([t.strip().upper() for t in targets])
    else:
        targets = list(ORFs.keys())

    return targets

def translate(seq, name) -> str:
    """Takes a nucleotide sequence and returns the translated protein sequence.

    Args:
        seq (str): This is the nucleotide sequence to translate
        name (str): The ORF name of the nucleotide sequence.

    Returns:
        str: The translated protein sequence.
    """
      
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if 'N' in codon:
                n_count = seq.count('N')
                print(f'DNA Sequence {name} contains {n_count} Ns! Skipping.')
                protein = ""
                return protein
            protein+= table[codon]
    else:
        print(f'DNA sequence {name} not a a multiple of 3!')
      
    return protein

def write_sbatch_file(seq_name, output_dir, AF_out) -> str:
    """Generates a BATCH file that will be used to launch the analysis on the Institut Pasteur HPC.
    As of March, 2022, the .sh file is compatible with the Alphafold module.

    Args:
        seq_name (str): The ORF name. Used to identify the protein, and the results.
        output_dir (str): The given output directory, to be prepended with the current absolute path.
        AF_out (str): The output directory that the batch file assigns to future Alphafold results.

    Returns:
        str: The formatted BATCH file.
    """
    cwd = os.getcwd()

    #  Unless we are doing some historical modeling, we use everything currently available (today)
    max_model_date = datetime.now().strftime('%Y-%m-%d')
    
    # Some additional control of where to put the eventual AF results was added.
    if AF_out:
        output_dir_str = f'OUTPUT_DIR={AF_out}/{seq_name}\n'
    else:
        output_dir_str = f'OUTPUT_DIR={cwd}/af_results/{seq_name}\n'

    text =  '#!/bin/bash\n' \
            '#SBATCH -c 12\n' \
            '#SBATCH -p common\n\n' + \
            f'#SBATCH -J {seq_name} #Job Name\n\n' \
            "OPENMM_PLATFORM='CPU'\n" \
            'OPENMM_CPU_THREADS=12\n' \
            'ALPHAFOLD_JACKHMMER_N_CPU=12\n' \
            'ALPHAFOLD_HHBLITS_N_CPU=12\n' \
            f'INPUT_FASTA={output_dir}/{seq_name}.aa.fasta\n' \
            f'{output_dir_str}' \
            'PRESET_MODEL=monomer\n' \
            'PRESET_DB=full_dbs\n' \
            'alphafold --fasta_paths ${INPUT_FASTA} \\\n' \
            '          --output_dir ${OUTPUT_DIR} \\\n' \
            f'          --max_template_date {max_model_date} \\\n' \
            '          --use_gpu_relax=false \\\n' \
            '          --model_preset ${PRESET_MODEL} \\\n' \
            '          --data_dir ${ALPHAFOLD_DATA} \\\n' \
            '          --uniref90_database_path ${ALPHAFOLD_DATA}/uniref90/uniref90.fasta \\\n' \
            '          --mgnify_database_path ${ALPHAFOLD_DATA}/mgnify/mgy_clusters.fa \\\n' \
            '          --pdb70_database_path ${ALPHAFOLD_DATA}/pdb70/pdb70 \\\n' \
            '          --template_mmcif_dir ${ALPHAFOLD_DATA}/pdb_mmcif/mmcif_files \\\n' \
            '          --obsolete_pdbs_path ${ALPHAFOLD_DATA}/pdb_mmcif/obsolete.dat \\\n' \
            '          --bfd_database_path ${ALPHAFOLD_DATA}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\\n' \
            '          --uniclust30_database_path ${ALPHAFOLD_DATA}/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \\\n' \
            '          --db_preset ${PRESET_DB} \\\n' \
            '          --verbosity 1'
    return text

def make_orfs(seqFile, format) -> pyfaidx.Fasta:
    """If necessary, depending on the format of the header of the FASTA seqFile, creates a reformatted temporary file
    that allows for simpler processing and selection of ORFs. 

    Args:
        seqFile (str): Multi FASTA file with all ORFs from the Genome
        format (str): The origin/format of the FASTA file ['refseq', 'genbank', 'none']

    Returns:
        pyfaidx.Fasta: Indexed list of ORFs with simplified names.
    """
    if format == 'none':
        return pyfaidx.Fasta(seqFile)
    elif format == 'refseq':
        orf_header = re.compile('^>\S+.(\S+.\[locus_tag=(\S+)].+)')
        tempfile = seqFile + '.tmp'
        t = open(tempfile, 'w+')

        with open(seqFile, 'r') as i:
            for line in i.readlines():
                if '>' in line:
                    match = orf_header.match(line)
                    t.write('>' + match.group(2) + ' ' + match.group(1) + '\n')
                else:
                    t.write(line)
        t.close()
        ORFs = pyfaidx.Fasta(tempfile)
        os.remove(tempfile)
        os.remove(tempfile + '.fai')
        return ORFs
    
    elif format == 'gff':
        seqFile_as_gff = gp.read_gff3(seqFile)
        print('Parsing sequence file in GFF format.\n%s' % seqFile_as_gff.header)
        df = seqFile_as_gff.df  # Just the entries, without header information

        df = df[(df['type'] == 'CDS')]  # We only want the coding sequences.
        print('%i ORFs found in GFF file.' % df.shape[0])

        tempfile = seqFile + '.tmp'
        t = open(tempfile, 'w+')        
        for i in df.index:  # For each CDS
            attributes = df.loc[i, 'attributes']  # attributes column, has info and sequence.
            attributes = attributes.split(';')  # Get a list of items ['attribute=value']
            
            for i in attributes:  # Inside each CDS's attributes
                if 'locus_tag' in i:  # Find the attribute that names the locus
                    t.write('>' + i.split('=')[1] + '\n')  # Get just the name
                if 'translation' in i: # Get the translation
                    t.write(i.split('=')[1]+ '\n')  # Get the translation
        t.close()

        ORFs = pyfaidx.Fasta(tempfile)
        os.remove(tempfile)
        os.remove(tempfile + '.fai')
        return ORFs

    else:
        pass

def is_dna(sequence_string) -> bool:
    dna = "ATCGN"
    sequence_string = sequence_string.upper()

    answer = all(i in dna for i in sequence_string)    
    return answer  # True=DNA, False=Protein

def run():
    args = parse_arguments()
    outdir = handle_outdir(args.outdir, args.skip_prompt)
    ORFs = make_orfs(args.seqFile, args.FASTA_Format)
    targets = make_targets(args.targetsFile, ORFs)

    for t in targets:
        nuc_seq = ORFs[t][:].seq
        if nuc_seq:  # does not run if sequence is empty
            if is_dna(nuc_seq):  # True = We confirm DNA and translate.
                aa_seq = f'>{t}\n' + translate(nuc_seq, t)[:-1].strip()
            else:  # False = We assume protein and continue
                aa_seq = nuc_seq        
        
            with open(outdir + '/' + t + '.aa.fasta', 'w+') as outfile:
                outfile.write(aa_seq)
            
            sbatch_text = write_sbatch_file(t, outdir, args.AF_out) 
            with open(outdir + '/' + t + '.sbatch', 'w+') as outfile:
                outfile.write(sbatch_text)

if __name__ == '__main__':
    run()
