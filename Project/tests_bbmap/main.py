import subprocess
import os
from Bio import SeqIO
import time
from pathlib import Path
import pandas as pd
from statistics import mode


def concat_references():
	sequences = []
	for path in Path(os.getcwd()).rglob('*.fasta'):
		if path.name != "primers.fasta" and path.name != "all_references.fasta":
			fasta_file = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
			file_name = path.name.replace('.fasta', '')
			
			seq_record = fasta_file[list(fasta_file.keys())[0]]
			seq_record.id = file_name
			seq_record.description = file_name
			sequences.append(seq_record)

	with open("all_references.fasta", "w") as output_handle:
		SeqIO.write(sequences, output_handle, "fasta")


def map2references(primer, primers):
	primer_code = primer.split('@')[0]
	primer_seq = primers[primer].seq

	output_file = os.path.join(os.getcwd(), primer_code, 'alignments.sam')

	# subprocess.run(
	# 	f'source /home/bioinfo/miniconda3/etc/profile.d/conda.sh && \
	# 	conda activate bbmap_env && \
	# 	msa.sh in=all_references.fasta out={output_file} literal={primer_seq} rcomp=t addr=t replicate=t cutoff=0.8 -Xmx30g',
	# 	shell=True,
	# 	executable='/bin/bash'
	# )

	alignments_df = pd.read_csv(output_file, sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])

	csv_output_file = os.path.join(os.getcwd(), primer_code, 'alignments.tsv')
	alignments_df.to_csv(csv_output_file, sep='\t', index=False)
	print(alignments_df.loc[alignments_df['RNAME'] != "*"])


def get_primers_alignment_data(primers_fasta):
	primers_alignment_data = {}
	for primer_name in primers_fasta.keys():
		primer_name = primer_name.split('@')[0]
		if primer_name not in primers_alignment_data.keys():
			primers_alignment_data[primer_name] = {'alignments': pd.read_csv(f'{primer_name}/alignments.tsv', sep='\t'), 'references': [file.replace('.fasta', '') for file in os.listdir(primer_name) if ".fasta" in file]}

	return primers_alignment_data


def get_common_organisms(alignments_0_df, alignments_1_df):		
	organisms_0 = set(alignments_0_df['RNAME'])
	
	organisms_1 = set(alignments_1_df['RNAME'])

	common_organisms = organisms_0.intersection(organisms_1)
	common_organisms.remove('*')

	return common_organisms


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
	elif pos_0 > pos_1 and query_0[0:2] != "r_" and query_1[0:2] != "r_":
		orientation = (pair[1], pair[0])
	elif pos_0 < pos_1 and query_0[0:2] == "r_" and query_1[0:2] != "r_":
		orientation = (pair[1], pair[0])
	elif pos_0 < pos_1 and query_0[0:2] == "r_" and query_1[0:2] == "r_":
		orientation = (pair[0], pair[1])
	elif pos_0 > pos_1 and query_0[0:2] == "r_" and query_1[0:2] == "r_":
		orientation = (pair[1], pair[0])

	return orientation


def get_primers_pairs(primers_alignment_data):
	primer_pairs = [(primer_a, primer_b) for idx, primer_a in enumerate(primers_alignment_data.keys()) for primer_b in list(primers_alignment_data.keys())[idx + 1:]]

	accepted_pairs = []
	for pair in primer_pairs:
		alignments_0_df = primers_alignment_data[pair[0]]['alignments']
		alignments_1_df = primers_alignment_data[pair[1]]['alignments']
		common_organisms = get_common_organisms(alignments_0_df, alignments_1_df)
		
		orientations = []
		if len(common_organisms) >= 10:  # pelo menos metade das referencias são comun entre eles 
			for organism in common_organisms:
				orientations.append(get_pair_orientation(alignments_0_df, alignments_1_df, pair, organism))
					
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
	start = time.process_time()

	# concat_references()

	primers_fasta = SeqIO.to_dict(SeqIO.parse("primers.fasta", "fasta"))

	# for primer in primers_fasta.keys():
	# 	map2references(primer, primers_fasta)


	primers_alignment_data = get_primers_alignment_data(primers_fasta)
	accepted_pairs = get_primers_pairs(primers_alignment_data)

	primers_meta = pd.read_csv('primers.csv', sep=';')
	write_pairs_file(accepted_pairs, primers_meta)
	

	end = time.process_time()
	total = end - start

	print('Tempo de execução: ', total)
			
			
# para corrigir problema com primers que se alinham em mais de uma localidade, posso fazer o alinhamento mais de uma vez, cortando a sequencia onde encontrar alinhamento, e terminar quando nao conseguir alinhar mais

