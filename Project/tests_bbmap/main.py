import subprocess
import os
from Bio import SeqIO
import time
from pathlib import Path
import pandas as pd
from statistics import mode
import concurrent.futures
import threading


def concat_references(file_size_cutoff):
	
	sequences = []
	references_file_size = 0
	references_files_count = 0

	directory = Path('primers_references')

	# Obtém todos os arquivos com a extensão .fasta no diretório e subdiretórios
	fasta_files = list(directory.rglob('*.fasta'))

	# Ordena os arquivos por tamanho
	# sorted_fasta_files = sorted(fasta_files, key=lambda file: file.stat().st_size)

	# Agora você pode iterar sobre os arquivos ordenados por tamanho
	for path in fasta_files:
		if path.name != "primers.fasta" and path.name != "all_references.fasta":
			fasta_file = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
			file_name = path.name.replace('.fasta', '')
			seq_record = fasta_file[list(fasta_file.keys())[0]]
			seq_record.id = str(path.parent).split('/')[-1] + '_' + file_name
			seq_record.description = ''

			sequences.append(seq_record)
			references_file_size += (1 + len(str(path.parent).split('/')[-1] + '_' + file_name) + 1 + len(seq_record.seq) + len(seq_record.seq) // 60 + 1) / 1000
			# if references_file_size / 1000 > file_size_cutoff:
			with open(f"all_references/all_references_{references_files_count}.fasta", "w") as output_handle:
				SeqIO.write(sequences, output_handle, "fasta")

			references_files_count += 1
			sequences = []

			#references_file_size = (1 + len(str(path.parent).split('/')[-1] + '_' + file_name) + 1 + len(seq_record.seq) + len(seq_record.seq) // 60 + 1) / 1000

			#else:
				#sequences.append(seq_record)
	# remaining_filename = 'all_references/all_references.fasta' if references_files_count == 0 else f"all_references/all_references_{references_files_count}.fasta"
	# with open(remaining_filename, "w") as output_handle:
	# 	SeqIO.write(sequences, output_handle, "fasta")


def map2references(primer, primers, reference_file, pos, memory):
	primer_seq = primers.loc[primers['id'] == primer]['nuc'].iloc[0]

	ambiguous = len(primer_seq) - primer_seq.count('A') - primer_seq.count('C') - primer_seq.count('T') - primer_seq.count('G') 
	if str(primer) in os.listdir('primers_references') and ambiguous <= 3:
	
		primer_seq = primers.loc[primers['id'] == primer]['nuc'].iloc[0]

		output_file = os.path.join('alignments', f'{primer}_{pos}.sam')


		# memory = int(0.3 * (os.path.getsize(f'all_references/{references_file}') / 1000000))
		subprocess.run(
			f'source /home/bioinfo/miniconda3/etc/profile.d/conda.sh && \
			conda activate bbmap_env && \
			msa.sh in=all_references/{references_file} out={output_file} literal={primer_seq} rcomp=t addr=t replicate=t cutoff=0.8 -Xmx{memory}g && \
			cat alignments/{primer}_*.sam >> alignments/{primer}.sam && \
			rm alignments/{primer}_*.sam',
			shell=True,
			executable='/bin/bash',
			stdout=subprocess.DEVNULL,
    		stderr=subprocess.STDOUT
		)


		# alignments_df = pd.read_csv(output_file, sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])

		# csv_output_file = os.path.join('alignments', f'{primer}_{pos}.tsv')
		# alignments_df.to_csv(csv_output_file, sep='\t', index=False)
	

def get_primers_alignment_data(primers_fasta):
	primers_alignment_data = {}
	for primer_id in primers_fasta.keys():
		primer_id = primer_id.split(' ')[0]
		if primer_id in os.listdir('primers_references'):
			if primer_id not in primers_alignment_data.keys():
				alignments_df = pd.read_csv(f'alignments/{primer_id}.sam', sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])
				# csv_output_file = os.path.join('alignments', f'{primer_id}.tsv')
				# alignments_df.to_csv(csv_output_file, sep='\t', index=False)
				# primer_alignments = [pd.read_csv(str(alignment), sep='\t') for alignment in Path('alignments').rglob(f'{primer_id}_*.tsv')]

				# alignments_df = pd.concat(primer_alignments).reset_index(drop=True)
				primers_alignment_data[primer_id] = {'alignments': alignments_df, 'references': [file.replace('.fasta', '') for file in os.listdir(f'primers_references/{primer_id}') if ".fasta" in file]}

	return primers_alignment_data


def get_common_organisms(alignments_0_df, alignments_1_df):		
	organisms_0 = set(alignments_0_df['RNAME'])
	
	organisms_1 = set(alignments_1_df['RNAME'])

	common_organisms = organisms_0.intersection(organisms_1)
	common_organisms.remove('*')

	return common_organisms


def check_product_size(alignments_0_df, alignments_1_df, organism):
	organism_0 = alignments_0_df.loc[alignments_0_df['RNAME'] == organism].reset_index(drop=True)
	pos_0 = organism_0['POS'][0]

	organism_1 = alignments_1_df.loc[alignments_1_df['RNAME'] == organism].reset_index(drop=True)
	pos_1 = organism_1['POS'][0]

	if abs(pos_0 - pos_1) >= 200 and abs(pos_0 - pos_1) <= 10000:
		return True
	else:
		return False


def get_pair_orientation(alignments_0_df, alignments_1_df, pair, organism):
	organism_0 = alignments_0_df.loc[alignments_0_df['RNAME'] == organism].reset_index(drop=True)
	pos_0 = organism_0['POS'][0]
	query_0 = organism_0['QNAME'][0]

	organism_1 = alignments_1_df.loc[alignments_1_df['RNAME'] == organism].reset_index(drop=True)
	pos_1 = organism_1['POS'][0]
	query_1 = organism_1['QNAME'][0]

	orientation = None
	if pos_0 < pos_1 and query_0[0:2] != "r_" and query_1[0:2] != "r_" :
		orientation = (pair[0], pair[1])
	elif pos_0 > pos_1 and query_0[0:2] != "r_" and query_1[0:2] == "r_":
		orientation = (pair[0], pair[1])
	elif pos_0 < pos_1 and query_0[0:2] != 'r_' and query_1[0:2] == 'r_':
		orientation = (pair[1], pair[0])
	elif pos_0 > pos_1 and query_0[0:2] != "r_" and query_1[0:2] != "r_":
		orientation = (pair[1], pair[0])
	elif pos_0 < pos_1 and query_0[0:2] == "r_" and query_1[0:2] != "r_":
		orientation = (pair[1], pair[0])
	elif pos_0 > pos_1 and query_0[0:2] == "r_" and query_1[0:2] != "r_":
		orientation = (pair[0], pair[1])
	elif pos_0 < pos_1 and query_0[0:2] == "r_" and query_1[0:2] == "r_":
		orientation = (pair[0], pair[1])
	elif pos_0 > pos_1 and query_0[0:2] == "r_" and query_1[0:2] == "r_":
		orientation = (pair[1], pair[0])

	else:
		print(pos_0, pos_1, query_0, query_1)
	return orientation


def get_primers_pairs(primers_alignment_data):
	primer_pairs = [(primer_a, primer_b) for idx, primer_a in enumerate(primers_alignment_data.keys()) for primer_b in list(primers_alignment_data.keys())[idx + 1:]]

	accepted_pairs = []
	for pair in primer_pairs:
		alignments_0_df = primers_alignment_data[pair[0]]['alignments']
		alignments_1_df = primers_alignment_data[pair[1]]['alignments']
		common_organisms = get_common_organisms(alignments_0_df, alignments_1_df)
		orientations = []
		if len(common_organisms) >= 20:  # pelo menos metade das referencias são comun entre eles 
			for organism in common_organisms:
				if check_product_size(alignments_0_df, alignments_1_df, organism):
					pair_orientation = get_pair_orientation(alignments_0_df, alignments_1_df, pair, organism)
					if pair_orientation is not None:
						orientations.append(pair_orientation)


		if len(orientations) >= 20:
			accepted_pairs.append(mode(orientations))

	return accepted_pairs


def write_pairs_file(accepted_pairs, primers_meta):

	pairs_dict = {'NAME': [], 'LEFT_PRIMER': [], 'RIGHT_PRIMER': []}
	for pair in accepted_pairs:
		name = f'{pair[0]}_{pair[1]}'
		LEFT_PRIMER = primers_meta.loc[primers_meta['code'] == pair[0]]['nuc'].iloc[0]  # ver que tem primers que eu renomeei no outro programa.... isso na verdade tem que ser uma limpeza inclusa no programa de extração dos primers do bold
		RIGHT_PRIMER = primers_meta.loc[primers_meta['code'] == pair[1]]['nuc'].iloc[0]

		pairs_dict['NAME'].append(name)
		pairs_dict['LEFT_PRIMER'].append(LEFT_PRIMER)
		pairs_dict['RIGHT_PRIMER'].append(RIGHT_PRIMER)

	pd.DataFrame(pairs_dict).to_csv('primers_pairs.tsv', sep='\t', index=False)


if __name__ == '__main__':
	"""
	start = time.process_time()

	concat_references(100)

	primers_meta = pd.read_csv('primers.csv', sep=';')
	primers_fasta = SeqIO.to_dict(SeqIO.parse('primers.fasta', "fasta"))

	for primer in primers_meta['id']:

		primer_seq = primers_meta.loc[primers_meta['id'] == primer]['nuc'].iloc[0]

		ambiguous = len(primer_seq) - primer_seq.count('A') - primer_seq.count('C') - primer_seq.count('T') - primer_seq.count('G') 
		if str(primer) in os.listdir('primers_references') and ambiguous <= 3:
			map2references(primer, primers_meta, '30g')


	primers_alignment_data = get_primers_alignment_data(primers_fasta)
	accepted_pairs = get_primers_pairs(primers_alignment_data)
	print(accepted_pairs)
	print(len(accepted_pairs))

	print('\n\n\n')

	# write_pairs_file(accepted_pairs, primers_meta)
	

	end = time.process_time()
	total = end - start

	print('Tempo de execução: ', total)
			
	"""


	start = time.process_time()

	# concat_references(25)


	primers_meta = pd.read_csv('primers.csv', sep=';')
	primers_fasta = SeqIO.to_dict(SeqIO.parse('primers.fasta', "fasta"))

	from itertools import repeat
	from math import ceil


	for pos, references_file in enumerate(os.listdir('all_references')):
		print(pos)
		# output_file = os.path.join('alignments', f'{primer}_{pos}.sam')
		memory = ceil(0.3 * (os.path.getsize(f'all_references/{references_file}') / 1000000))

		if memory == 0:
			memory = 1
		threads = int(60/memory)

		if threads > 12:
			threads = 12

		if memory * threads < 60:
			memory = 60 / threads

		with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
			results = list(executor.map(map2references, primers_meta['id'], repeat(primers_meta), repeat(references_file), repeat(pos), repeat(memory)))

		

	# primers_alignment_data = get_primers_alignment_data(primers_fasta)
	# accepted_pairs = get_primers_pairs(primers_alignment_data)
	# print(accepted_pairs)
	# print(len(accepted_pairs))

	# print('\n\n\n')


		# write_pairs_file(accepted_pairs, primers_meta)


	end = time.process_time()
	total = end - start

	print('Tempo de execução: ', total)



# para corrigir problema com primers que se alinham em mais de uma localidade, posso fazer o alinhamento mais de uma vez, cortando a sequencia onde encontrar alinhamento, e terminar quando nao conseguir alinhar mais

# ver qual estrategia é melhor: usar o msa ou criar um banco blast
