import subprocess
import pandas as pd
from statistics import mode


def run_blast(output, max_target_seqs='12000', num_threads='10', database='/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/primers_references_nt'):
	# subprocess.run(
	# 	f'export BLASTDB=/media/bioinfo/6tb_hdd/04_Blast_Databases/taxdb && \
	# 	time blastn \
	# 	-db {database} \
	# 	-query /media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/primers_seqs.fa \
	# 	-outfmt "6 qseqid pident gaps qcovs mismatch evalue score bitscore sgi sacc slen qstart qend sstart send qseq sseq length sstrand frames stitle sscinames" \
	# 	-out {output} \
	# 	-num_threads 10 \
	# 	-max_target_seqs 500 \
	# 	-task blastn-short',
	# 	shell=True,
	# 	executable='/bin/bash'
	# )

	subprocess.run(
		f'export BLASTDB=/media/bioinfo/6tb_hdd/04_Blast_Databases/taxdb && \
		time blastn \
		-db {database} \
		-query /media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/primers_seqs.fa \
		-outfmt "6 qseqid frames sstart sacc send length mismatch" \
		-out {output} \
		-num_threads 10 \
		-max_target_seqs 500 \
		-task blastn-short',
		shell=True,
		executable='/bin/bash'
	)

	


def split_and_duplicate_primer(row, sep):
    if sep in row:
        splitted_values = row.split(sep)
        primer = ''.join(splitted_values[0:-2])
        pos = splitted_values[-1]
        return pd.Series([splitted_values[0], pos])
    else:
        return pd.Series([row, row])


def split_orientation(row, sep):
    if sep in row:
        splitted_values = row.split(sep)
        return splitted_values[-1]


def get_primers_alignment_data(primers, df):
	primers_alignment_data = {}
	for primer_id in primers:

		primer_df = df.loc[df['Primer'] == primer_id].reset_index(drop=True)

		primers_alignment_data[primer_id] = primer_df

	return primers_alignment_data


def clean_primers_alignment_data(primers_alignment_data, primers):
	for item in primers_alignment_data.items():
		primer_id = item[0]
		primer_df = item[1]
		primer_seq = primers.loc[primers['id'] == primer_id]['nuc'].iloc[0]
		
		# alinhamentos do tamanho do primer até o tamanho do primer - 2
		# primer_df = primer_df.loc[primer_df['length'].isin((len(primer_seq), len(primer_seq) - 1, len(primer_seq) - 2))]
		primer_df = primer_df.loc[primer_df['length'] == len(primer_seq)]

		# até 2 mismatches
		# primer_df = primer_df.loc[primer_df['mismatch'] <= 2]
		primer_df = primer_df.loc[primer_df['mismatch'] == 0]

		primers_alignment_data[primer_id] = primer_df

	return primers_alignment_data


def get_common_organisms(alignments_0_df, alignments_1_df):		
	organisms_0 = set(alignments_0_df['sacc'])
	
	organisms_1 = set(alignments_1_df['sacc'])

	common_organisms = organisms_0.intersection(organisms_1)
	return common_organisms


def check_product_size(alignments_0_df, alignments_1_df, organism):
	organism_0 = alignments_0_df.loc[alignments_0_df['sacc'] == organism].reset_index(drop=True)
	pos_0 = organism_0['sstart'][0]

	organism_1 = alignments_1_df.loc[alignments_1_df['sacc'] == organism].reset_index(drop=True)
	pos_1 = organism_1['send'][0]

	if abs(pos_0 - pos_1) >= 500 and abs(pos_0 - pos_1) <= 3000:
		return True
	else:
		return False


def get_pair_orientation(alignments_0_df, alignments_1_df, pair, organism, primers):
	direction_0 = primers.loc[primers['id'] == pair[0]]['direction'].iloc[0]
	direction_1 = primers.loc[primers['id'] == pair[1]]['direction'].iloc[0]

	organism_0 = alignments_0_df.loc[alignments_0_df['sacc'] == organism].reset_index(drop=True)
	pos_0 = organism_0['sstart'][0]
	orientation_0 = organism_0['orientation'][0]

	organism_1 = alignments_1_df.loc[alignments_1_df['sacc'] == organism].reset_index(drop=True)
	pos_1 = organism_1['sstart'][0]
	orientation_1 = organism_1['orientation'][0]

	smallest_position = min((pos_0, pos_1)) 

	orientation = None

	# orientation = (pair[0], pair[1]) if pos_0 < pos_1 else (pair[1], pair[0])
	if direction_0 != direction_1:
		if direction_0 == 'F':
			orientation = (pair[0], pair[1])
		else:
			orientation = (pair[1], pair[0])
	else:
		if pos_1 < pos_0 and orientation_1 == '-1' and orientation_0 == '1':
			orientation = (pair[1], pair[0])
		elif pos_0 < pos_1 and orientation_0 == '-1' and orientation_1 == '1':
			orientation = (pair[0], pair[1])
		elif pos_1 < pos_0 and orientation_1 == '1' and orientation_0 == '-1':
			orientation = (pair[1], pair[0])
		elif pos_0 < pos_1 and orientation_0 == '1' and orientation_1 == '-1':
			orientation = (pair[0], pair[1])

	return orientation


def get_primers_pairs(primers_alignment_data, primers):
	primer_pairs = [(primer_a, primer_b) for idx, primer_a in enumerate(primers_alignment_data.keys()) for primer_b in list(primers_alignment_data.keys())[idx + 1:]]

	print(f'Total de pares a serem processados: {len(primer_pairs)}')
	print('\n\n\n')
	accepted_pairs = []
	for pos, pair in enumerate(primer_pairs):
		common_organisms = get_common_organisms(primers_alignment_data[pair[0]], primers_alignment_data[pair[1]])

		orientations = []
		if len(common_organisms) >= 15:  # pelo menos metade das referencias são comun entre eles 
			for organism in common_organisms:
				check = check_product_size(primers_alignment_data[pair[0]], primers_alignment_data[pair[1]], organism)
				if check:
					pair_orientation = get_pair_orientation(primers_alignment_data[pair[0]], primers_alignment_data[pair[1]], pair, organism, primers)
					if pair_orientation is not None:
						orientations.append(pair_orientation)

		if len(orientations) >= 15:
			print(f'\tPar {pos} aceito, orientação:', mode(orientations))
			accepted_pairs.append(mode(orientations))

	return accepted_pairs


def write_pairs_file(accepted_pairs, primers_meta):
	pairs_dict = {'NAME': [], 'LEFT_PRIMER': [], 'RIGHT_PRIMER': []}
	for pair in accepted_pairs:
		name = f'{pair[0]}_{pair[1]}'
		LEFT_PRIMER = primers_meta.loc[primers_meta['id'] == pair[0]]['nuc'].iloc[0]  # ver que tem primers que eu renomeei no outro programa.... isso na verdade tem que ser uma limpeza inclusa no programa de extração dos primers do bold
		RIGHT_PRIMER = primers_meta.loc[primers_meta['id'] == pair[1]]['nuc'].iloc[0]

		pairs_dict['NAME'].append(name)
		pairs_dict['LEFT_PRIMER'].append(LEFT_PRIMER)
		pairs_dict['RIGHT_PRIMER'].append(RIGHT_PRIMER)

	pd.DataFrame(pairs_dict).to_csv('500_customdb_less_tolerance_15organisms.tsv', sep='\t', index=False)


if __name__ == '__main__':

	# run_blast('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_custom_db_500.tsv', database='/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test1/db/references')
	# run_blast('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_nt_500.tsv')
	# run_blast('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_16s_500.tsv', database='/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/primers_references_16s')

	
	primers = pd.read_csv('primers_filtered.tsv', sep='\t', dtype={'id': str})

	print('Lendo arquivo BLAST')
	# df = pd.read_csv('/media/bioinfo/6tb_hdd/03_ELLEN/02_data/primers_references_blast_alignment/blast_output.tsv', sep='\t', names=['qseqid' ,'pident' ,'gaps' ,'qcovs' ,'mismatch' ,'evalue' ,'score' ,'bitscore' ,'sgi' ,'sacc' ,'slen' ,'qstart' ,'qend' ,'sstart' ,'send' ,'qseq' ,'sseq' ,'length' ,'sstrand' ,'frames' ,'stitle' ,'sscinames'], dtype={'qseqid': str})

	df_nt = pd.read_csv('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_nt.tsv', sep='\t', names=['qseqid' ,'pident' ,'gaps' ,'qcovs' ,'mismatch' ,'evalue' ,'score' ,'bitscore' ,'sgi' ,'sacc' ,'slen' ,'qstart' ,'qend' ,'sstart' ,'send' ,'qseq' ,'sseq' ,'length' ,'sstrand' ,'frames' ,'stitle' ,'sscinames'], dtype={'qseqid': str})
	df_nt = df_nt[['qseqid', 'frames', 'sstart', 'sacc', 'send', 'length', 'mismatch']]
	# df_nt = pd.read_csv('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_nt.tsv', sep='\t', names=['qseqid', 'frames', 'sstart', 'sacc', 'send', 'length', 'mismatch'], dtype={'qseqid': str})

	df_16s = pd.read_csv('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_16s.tsv', sep='\t', names=['qseqid' ,'pident' ,'gaps' ,'qcovs' ,'mismatch' ,'evalue' ,'score' ,'bitscore' ,'sgi' ,'sacc' ,'slen' ,'qstart' ,'qend' ,'sstart' ,'send' ,'qseq' ,'sseq' ,'length' ,'sstrand' ,'frames' ,'stitle' ,'sscinames'], dtype={'qseqid': str})

	df_16s = df_16s[['qseqid', 'frames', 'sstart', 'sacc', 'send', 'length', 'mismatch']]
	# df_16s = pd.read_csv('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_16s.tsv', sep='\t', names=['qseqid', 'frames', 'sstart', 'sacc', 'send', 'length', 'mismatch'], dtype={'qseqid': str})


	df = pd.concat([df_nt, df_16s]).reset_index(drop=True)

	df = pd.read_csv('/media/bioinfo/6tb_hdd/03_ELLEN/test_blastdb_gilist/test2/blast_output_custom_db_500.tsv', sep='\t', names=['qseqid', 'frames', 'sstart', 'sacc', 'send', 'length', 'mismatch'], dtype={'qseqid': str})
	# df_test = df.head(200000)
	# df_test.to_csv('blast_test.tsv', sep='\t', index=False)


	# df = pd.read_csv('blast_test.tsv', sep='\t', dtype={'qseqid': str})
	# df[['Primer', 'SeqID']] = df['qseqid'].apply(split_and_duplicate_primer, args='_')

	print('Separando primers e posições')
	df[['Primer', 'SeqID']] = df['qseqid'].str.split('_', n=1, expand=True)
	print('Separando orientações')
	df[['tmp', 'orientation']] = df['frames'].str.split('/', n=1, expand=True)

	print('Obtendo dataframes dos primers')
	primers_alignment_data = get_primers_alignment_data(df['Primer'].unique(), df)
	print('Limpando dataframes dos primers')
	primers_alignment_data = clean_primers_alignment_data(primers_alignment_data, primers)
	
	print('Obtendo pares de primers')
	accepted_pairs = get_primers_pairs(primers_alignment_data, primers)
	print('Salvando pares de primers')
	write_pairs_file(accepted_pairs, primers)
