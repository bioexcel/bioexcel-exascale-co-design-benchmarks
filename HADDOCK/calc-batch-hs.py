import glob
import sys
import re

"""
 Calculates HADDOCK-score for a list of PDB structures
 
Author: Rodrigo Honorato
Date: 10-July-2019
"""


def create_targetlist(path):
	""" Get the location of all PDB decoys  """
	target_list = glob.glob(path + '/*.pdb')
	return target_list


def read_targetfile(f):
	""" Read file list and return a list """
	target_list = []
	for e in open(f).readlines():
		pdb = e.rstrip()
		target_list.append(pdb)
	return target_list


def extract_energies(pdb_file):
	""" Extract energies from the header of the PDB file, according to HADDOCK formatting """
	vdw = .0
	elec = .0
	desolv = .0
	air = .0
	bsa = .0

	vdw_elec_air_regex = r"\s(\-?\d*\.?\d{1,}|0\b)"  # Known error: 8.754077E-02 is matched as 8.754077 which is
	# not relevant for the I/O benchmark but should be taken into account

	desolv_regex = r"(\-?\d*\.?\d*)$"
	bsa_regex = r"(\-?\d*\.?\d*)$"

	f = open(pdb_file, 'r')
	for line in f:

		if 'REMARK energies' in line:
			# print(line)
			total, bonds, angles, improper, dihe, vdw, elec, air, cdih, coup, rdcs, vean, dani, xpcs, rg = re.findall(
				vdw_elec_air_regex, line)
			vdw = float(vdw)
			elec = float(elec)
			air = float(elec)

		if 'REMARK Desolvation' in line:
			# print(line)
			desolv = float(re.findall(desolv_regex, line)[0])

		if 'REMARK buried surface area' in line:
			# print(line)
			bsa = float(re.findall(bsa_regex, line)[0])
			break

	f.close()

	return vdw, elec, desolv, air, bsa


def calculate_haddock_score(pdb_file):
	""" Calculate the HADDOCK Score of the PDB file using its appropriate weight """

	weight_dic = dict(it0=(0.01, 1.0, 1.0, 0.01, 0.01), it1=(1.0, 1.0, 1.0, 0.1, 0.01), itw=(1.0, 0.2, 1.0, 0.1, 0.0))

	vdw, elec, desolv, air, bsa = extract_energies(pdb_file)
	phase_regex = r"-(it.)_"
	try:
		stage = re.findall(phase_regex, pdb_file)[0]
	except IndexError:
		stage = 'it0'

	w = weight_dic[stage]
	haddock_score = w[0] * vdw + w[1] * elec + w[2] * desolv + w[3] * air - w[4] * bsa

	return haddock_score


def output(result_list):
	result_list.sort(key=lambda x: x[1])
	with open('file.list', 'w') as f:
		for decoy_name, haddock_score in result_list:
			f.write(f'{decoy_name} {haddock_score:.4f}\n')
	f.close()


def main():
	if len(sys.argv) > 1:
		target_f = sys.argv[1]
	else:
		print('Please specify the file list with the PDB paths')
		sys.exit(1)

	file_list = []
	pdb_list = read_targetfile(target_f)
	# pdb_list = create_targetlist(target_f)

	for decoy in pdb_list:
		haddock_score = calculate_haddock_score(decoy)
		file_list.append((decoy, haddock_score))

	output(file_list)


if __name__ == '__main__':
	main()
