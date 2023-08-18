import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Align


primers = pd.read_csv('primers.csv', sep=';')

def count_ambiguous(seq):
	iupac_table = {
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
	count = 0
	for base in seq:
		if base in iupac_table.keys():
			count += 1

	return count

primers_sequences = []
for index, primer in primers.iterrows():
	ambiguous_count = count_ambiguous(primer['nuc'])
	if ambiguous_count <= 2:
		new_record = SeqRecord(
			Seq(primer['nuc']).reverse_complement(),
			id = str(primer['id']),
			description='primer'
		)

		primers_sequences.append(new_record)

SeqIO.write(primers_sequences, 'primers.fasta', "fasta")

primers_dict = {}
for primer in primers_sequences:
	primers_dict[primer.id] = primer.seq


def get_mismatches(reference_sequence, matched_sequence, pos, primer):
	iupac_table = {
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


	alignment = reference_sequence.seq[pos: pos + len(matched_sequence)]


	if len(matched_sequence) < len(primer):
		# print('menor')
		for i in range(len(primer) - len(matched_sequence)):
			matched_sequence += '-'
			alignment += '-'


	mismatches = len(primer)
	for i, base1 in enumerate(matched_sequence):
	    base2 = alignment[i]
	    if 'I' in base1 or 'I' in base2:
	    	mismatches = '?I'
	    	break
	    elif '-' not in base2:
		    for iupac_base in iupac_table[base1]:
		    	for iupac_base2 in iupac_table[base2]:
				    if iupac_base == iupac_base2:
				        mismatches -= 1

	# print(f'{mismatches}/{len(primer)}')
	return mismatches


results_dict = {'Primer': []}
repeated_primers_mismatches = {'Primer': []}

count = 0
for reference in os.listdir('references'):
	if ' ' in reference:
		reference = reference.replace(' ', '_')
		os.rename(f'references/{reference}', f'references/{reference}')

	basename = reference.replace('.fasta', '')
	basename = basename.replace('.fa', '')
	basename = basename.replace('.fna', '')
	print(basename)
	subprocess.run(
		f'mafft --anysymbol --keeplength --addfragments primers.fasta references/{reference} > results/{basename}_alignment.fasta',
		shell=True,
		executable='/bin/bash'
	)

	alignment_file = f'results/{basename}_alignment.fasta'

	alignment = Align.MultipleSeqAlignment(SeqIO.parse(alignment_file, "fasta"))
	reference_sequence = alignment[0]
	aligned_sequences = alignment[1:]

	for sequence in aligned_sequences:
		alignment_start_pos = -1

		for i, base in enumerate(sequence.seq):
			if base != "-":
				alignment_start_pos = i + 1

				if sequence.id not in results_dict['Primer']:
					results_dict['Primer'].append(sequence.id)

				if basename not in list(results_dict.keys()):
					results_dict[basename] = []

				
				break

		if alignment_start_pos != -1:
			matched_sequence = ''
			for base in sequence[i:]:
				if base != "-":
					matched_sequence += base
				else:
					break
			mismatches = get_mismatches(reference_sequence, matched_sequence, i, primers_dict[sequence.id])

			if mismatches != '?I' and mismatches <= 1:
				results_dict[basename].append(alignment_start_pos)
			else:
				results_dict[basename].append(None)

		if alignment_start_pos == -1:
			if sequence.id not in results_dict['Primer']:
				results_dict['Primer'].append(sequence.id)
			if basename not in list(results_dict.keys()):
				results_dict[basename] = []
			results_dict[basename].append(None)

	
print(repeated_primers_mismatches)
df = pd.DataFrame(results_dict)
df.to_csv('results/results6_revcomp_less_degenerated.tsv', sep='\t', index=False)
