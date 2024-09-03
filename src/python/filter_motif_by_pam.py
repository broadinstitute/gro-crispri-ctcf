'''
Read the MEME motifs, extract the PFM matrix and select the motifs
that contain a PAM sequence.

argument 1 is the folder with the meme motifs
argument 2 is the output file with the positive matches
argument 3 is the output file with the negative matches

#The meme motifs zip file (argument 1) can be downloaded from https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.zip

'''

from __future__ import division
import sys
from sys import argv
import glob


def read_meme_file(meme_file):
	'''
	Read the meme file and return the TF name and the 
	PFM matrix
	'''
	n = 0
	matrix = []
	for line in open(meme_file).read().split('\n')[:-1]:
		if line.startswith('MOTIF'):
			TF = line.split()[2]
		if line.startswith('URL'):
			n = 0
		if n == 1:
			line_matrix = line.split()
			matrix.append(line_matrix)
		if line.startswith('letter-probability'):
			n = 1
	return TF, matrix


def extract_patterns(matrix,  thr = 0.8):
	'''
	Find the instances in which there are two CCs in the first
	half of the motif or two GGs in the second half of the motif
	and report their location
	'''
	is_tf = []
	line_cnt  = 0
	number_pos = []; final_number_pos = []
	for line in matrix:
		if 0 < line_cnt < (len(matrix)/2+1):
			if (float(line[1])>thr):
				number_pos.append(line_cnt)
			else:
				number_pos = []
			if len(number_pos) > 1:
				is_tf.append('C')
				final_number_pos.append(number_pos)
		line_cnt += 1
	
	line_cnt  = 0
	number_pos = []
	for line in matrix:
		if (len(matrix)/2-1) < line_cnt < len(matrix):
			if (float(line[2])>thr):
				number_pos.append(line_cnt)
			else:
				number_pos = []
			if len(number_pos) > 1:
				is_tf.append('G')
				final_number_pos.append(number_pos)
		line_cnt += 1
	return is_tf, final_number_pos


def strong_motif(matrix, nt_list, number_pos_list, thr = 0.5):
	'''
	Extract the substring delimited by the PAM sequence and 
	select the motifs that have at least 5 nucleotides with 
	more than 0.5 probability in the PFM
	'''
	is_strong = False
	for nt in nt_list:
		pattern_strength = []
		if nt == 'C':
			number_pos = number_pos_list[nt_list.index('C')]
			matrix_mod = matrix[number_pos[0]:]
		if nt == 'G':
			number_pos = number_pos_list[nt_list.index('G')]
			matrix_mod = matrix[:number_pos[-1]+1]
		for line in matrix_mod:
			pattern_strength.append(max(line))
		counter = 0
		strong_n = 0; regular_n = 0

		if nt == 'G':
			for elem in list(reversed(pattern_strength)):
				if counter < 3:
					if float(elem) > 0.7:
						strong_n += 1
				else:
					if float(elem) > thr:
						regular_n += 1
				counter += 1

		if nt == 'C':
			for elem in pattern_strength:
				if counter < 3:
					if float(elem) > 0.7:
						strong_n += 1
				else:
					if float(elem) > thr:
						regular_n += 1
				counter += 1

		if (strong_n==3) & ((regular_n+strong_n)>5):
			is_strong = True

	return is_strong


if __name__ == "__main__":

	positives_file = open(argv[2],'w')
	negatives_file = open(argv[3],'w')

	for meme in glob.glob(argv[1]+"*meme"):

		TF, matrix = read_meme_file(meme)

		is_tf, number_pos = extract_patterns(matrix)

		if len(is_tf)>0:
			is_strong = strong_motif(matrix, is_tf, number_pos)

			if is_strong:
				positives_file.write(TF)
				positives_file.write('\n')

			if not is_strong:
				negatives_file.write(TF)
				negatives_file.write('\n')

		else:
			negatives_file.write(TF)
			negatives_file.write('\n')

	negatives_file.close()
	positives_file.close()






