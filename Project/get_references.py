import pandas as pd
import subprocess
from itertools import product
import os
import random
import time
import re
import yaml


def read_yaml(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


def extend_ambiguous_dna(seq):
	degenerated_table = {
		'A': ['A'],
		'C': ['C'],
		'G': ['G'],
		'T': ['T'],
		'W': ['A', 'T'],
		'S': ['C', 'G'],
		'M': ['A', 'C'],
		'K': ['G', 'T'],
		'R': ['A', 'G'],
		'Y': ['C', 'T'],
		'B': ['C', 'G', 'T'],
		'D': ['A', 'G', 'T'],
		'H': ['A', 'C', 'T'],
		'V': ['A', 'C', 'G'],
		'N': ['A', 'C', 'G', 'T'],
	}

	new_seq = ''
	for base in seq:
		if base == 'I':
			new_seq += 'N'
		else:
			new_seq += base
	seq = new_seq

	r = []
	# Aplica produto cartesiano nos conjuntos (conjunto = possíveis bases para cada letra)
	for i in product(*[degenerated_table[j] for j in seq]):
		r.append("".join(i))

	return r


def create_fasta(primers, output_file):
	seqs = []
	for primer_id, primer_code, primer_seq in zip(primers['id'], primers['code'], primers['nuc']):
		primer_seq_extended = extend_ambiguous_dna(primer_seq)
		try:
			if len(primer_seq_extended) == 1:
				seqs.append((str(primer_id), primer_seq_extended[0]))
			else:
				for pos, extended_primer in enumerate(primer_seq_extended):
					seqs.append((f'{primer_id}_{pos}', extended_primer))

		except:
			print(primer_code)

	with open(output_file, "w") as fasta_file:
		for primer_code, seq in seqs:
			fasta_record = f">{primer_code}\n{seq}\n"
			fasta_file.write(fasta_record)


def run_blast(fasta_file, output_file, database, tax_db, max_target_seqs='10', num_threads='1'):
	subprocess.run(
		f'export BLASTDB={tax_db} && \
		time blastn \
		-db {database} \
		-query {fasta_file} \
		-outfmt "6 qseqid pident gaps qcovs mismatch evalue score bitscore sgi sacc slen stitle sscinames" \
		-out {output_file} \
		-num_threads {num_threads} \
		-max_target_seqs {max_target_seqs} \
		-task blastn-short',
		shell=True,
		executable='/bin/bash'
	)


def split_and_duplicate(row):
    if '_' in row:
        valores_divididos = row.split('_')
        primer = ''.join(valores_divididos[0:-2])
        pos = valores_divididos[-1]
        return pd.Series([valores_divididos[0], pos])
    else:
        return pd.Series([row, row])


def clean_blast_results(df_references):
	df_references['sscinames'] = df_references['sscinames'].str.replace(' ', '_')
	df_references['sscinames'] = df_references['sscinames'].str.replace('(', '-')
	df_references['sscinames'] = df_references['sscinames'].str.replace(')', '-')

	full_pattern = re.compile('[^a-zA-Z0-9_-]')
	df_references['sscinames'] = df_references['sscinames'].str.replace(full_pattern, '_', regex=True)


	df_references[['Primer', 'SeqID']] = df_references['qseqid'].apply(split_and_duplicate)

	return df_references


def get_gis(df_references, primer_code, primer_seq, max_len, count, selected_organisms, path2save_references):
	df_primer_references = df_references.loc[df_references['Primer'] == primer_code].reset_index(drop=True)

	allowed_perc = (len(primer_seq) - 2) / len(primer_seq) * 100

	if max_len is not None:
		df_primer_references = df_primer_references.loc[df_primer_references['slen'] <= max_len].reset_index(drop=True)
	df_primer_references = df_primer_references.loc[df_primer_references['pident'] >= allowed_perc].reset_index(drop=True)
	df_primer_references = df_primer_references.loc[df_primer_references['qcovs'] >= allowed_perc].reset_index(drop=True)

	df_primer_references = df_primer_references.drop_duplicates(subset=['sscinames'])
	df_primer_references = df_primer_references.dropna(subset=['sscinames'])

	gis_nt = []
	gis_16S = []
	for organism, sgi, db in zip(df_primer_references['sscinames'], df_primer_references['sgi'], df_primer_references['level_0']):
		organism_renamed = organism.replace(' ', '_')
		try:
			if organism not in selected_organisms:
				selected_organisms.append(organism)
				count += 1
				if db == 'nt':
					gis_nt.append(str(sgi))
				elif db == '16S':
					gis_16S.append(str(sgi))

		except Exception as e:
			print(e)

		if count == 10:
			return count, selected_organisms, gis_nt, gis_16S

	return count, selected_organisms, gis_nt, gis_16S


def get_references(config, df_references, primers):
	df_references = df_references.sample(frac = 1).reset_index(drop=True)
	df_references = df_references.sort_values(by=['qcovs', 'pident', 'slen'], ascending=[False, False, False]).reset_index(drop=True)
	primer_codes = df_references['Primer'].unique()
	gis_nt_list = []
	gis_16s_list = []
	for primer_code in primer_codes:		
		count = 0
		max_len = (30000, 40000, 50000, 100000, 200000, None)
		len_pos = 0
		selected_organisms = []
		primer_seq = primers.loc[primers['id'] == primer_code]['nuc'].iloc[0]
		current_gis_nt = []
		current_gis_16S = []
		while count < 10 and len_pos < 6:
			count, selected_organisms, new_gis_nt, new_gis_16s = get_gis(df_references, primer_code, primer_seq, max_len[len_pos], count, selected_organisms, config['OUTPUT_PATH'])
			current_gis_nt.extend(new_gis_nt)
			current_gis_16S.extend(new_gis_16s)
			len_pos += 1


		gis_nt_list.extend(current_gis_nt)
		gis_16s_list.extend(current_gis_16S)

	return gis_nt_list, gis_16s_list


def create_alias(config, gis_nt_list, gis_16s_list):
	gis_nt = '\n'.join(gis_nt_list)
	with open(f"{config['OUTPUT_PATH']}gi_list_nt.txt", 'w') as f:
		f.write(gis_nt)

	gis_16S = '\n'.join(gis_16s_list)
	with open(f"{config['OUTPUT_PATH']}gi_list_16S.txt", 'w') as f:
		f.write(gis_16S)

	# blastdb_aliastool -gilist gi_list_nt.txt -db /media/bioinfo/6tb_hdd/04_Blast_Databases/new_nt/nt -out primers_references_nt -title "primers_references_nt" -dbtype nucl

	subprocess.run(
		f"cd {config['OUTPUT_PATH']} && \
		blastdb_aliastool -gilist gi_list_16S.txt -db {config['16S_DB_PATH']} -out primers_references_16s -title primers_references_16s -dbtype nucl",
		shell=True,
		executable='/bin/bash'
	)

	subprocess.run(
		f"cd {config['OUTPUT_PATH']} && \
		blastdb_aliastool -gilist gi_list_nt.txt -db {config['NT_DB_PATH']} -dbtype nucl -out primers_references_nt -title primers_references_nt",
		shell=True,
		executable='/bin/bash'
	)


def run():
	config = read_yaml('config.yaml')
	primers = pd.read_csv(config['PRIMERS'], sep='\t', dtype={'id': 'str'})
		
	create_fasta(primers, f"{config['OUTPUT_PATH']}primers_seqs.fa")

	# print('Rodando BLAST NT')
	# run_blast(f"{config['OUTPUT_PATH']}primers_seqs.fa", f"{config['OUTPUT_PATH']}primers_seqs_nt.tsv", database=config['NT_DB_PATH'], tax_db=config['TAX_DB_PATH'], max_target_seqs='50', num_threads=config['THREADS'])
	# print('Rodando BLAST 16S')
	# run_blast(f"{config['OUTPUT_PATH']}primers_seqs.fa", f"{config['OUTPUT_PATH']}primers_seqs_16S.tsv", database=config['16S_DB_PATH'], tax_db=config['TAX_DB_PATH'], max_target_seqs='50', num_threads=config['THREADS'])

	df_references_nt = pd.read_csv(f"{config['OUTPUT_PATH']}primers_seqs_nt.tsv", sep='\t', names=['qseqid', 'pident', 'gaps', 'qcovs', 'mismatch', 'evalue', 'score', 'bitscore', 'sgi', 'sacc', 'slen', 'stitle', 'sscinames'], dtype={'qseqid': 'str'})
	df_references_16S = pd.read_csv(f"{config['OUTPUT_PATH']}primers_seqs_16S.tsv", sep='\t', names=['qseqid', 'pident', 'gaps', 'qcovs', 'mismatch', 'evalue', 'score', 'bitscore', 'sgi', 'sacc', 'slen', 'stitle', 'sscinames'], dtype={'qseqid': 'str'})

	df_references = pd.concat([df_references_nt, df_references_16S], keys=['nt', '16S']).reset_index()
	df_references = clean_blast_results(df_references)

	gis_nt_list, gis_16s_list = get_references(config, df_references, primers)

	create_alias(config, gis_nt_list, gis_16s_list)



if __name__ == '__main__':
	config = read_yaml('config.yaml')
	primers = pd.read_csv(config['PRIMERS'], sep='\t', dtype={'id': 'str'})
		
	create_fasta(primers, f"{config['OUTPUT_PATH']}primers_seqs.fa")

	print('Rodando BLAST NT')
	run_blast(f"{config['OUTPUT_PATH']}primers_seqs.fa", f"{config['OUTPUT_PATH']}primers_seqs_nt.tsv", database=config['NT_DB_PATH'], tax_db=config['TAX_DB_PATH'], max_target_seqs='50', num_threads=config['THREADS'])
	print('Rodando BLAST 16S')
	run_blast(f"{config['OUTPUT_PATH']}primers_seqs.fa", f"{config['OUTPUT_PATH']}primers_seqs_16S.tsv", database=config['16S_DB_PATH'], tax_db=config['TAX_DB_PATH'], max_target_seqs='50', num_threads=config['THREADS'])

	df_references_nt = pd.read_csv(f"{config['OUTPUT_PATH']}primers_seqs_nt.tsv", sep='\t', names=['qseqid', 'pident', 'gaps', 'qcovs', 'mismatch', 'evalue', 'score', 'bitscore', 'sgi', 'sacc', 'slen', 'stitle', 'sscinames'], dtype={'qseqid': 'str'})
	df_references_16S = pd.read_csv(f"{config['OUTPUT_PATH']}primers_seqs_16S.tsv", sep='\t', names=['qseqid', 'pident', 'gaps', 'qcovs', 'mismatch', 'evalue', 'score', 'bitscore', 'sgi', 'sacc', 'slen', 'stitle', 'sscinames'], dtype={'qseqid': 'str'})

	df_references = pd.concat([df_references_nt, df_references_16S], keys=['nt', '16S']).reset_index()
	df_references = clean_blast_results(df_references)

	gis_nt_list, gis_16s_list = get_references(config, df_references, primers)

	create_alias(config, gis_nt_list, gis_16s_list)
