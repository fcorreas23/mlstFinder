import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord







def main():

	try:
		sys.argv[1]

	except Exception as e:
		print("ERROR: Faltan Parametros")
		sys.exit(0)

	

	print(f'Buscando MLST (acnB, gap, glta, gyrB, pgi, rpoD ) en { os.path.basename(sys.argv[1]) }')
	print('Ejecutando nhmmer....')
	
	subprocess.call(['nhmmer','--cpu', '6', '-E', '0.05' ,'-o', 'result_nhmmer.txt','--tblout', 'tmp.txt', 'hmm/mlts.hmm', sys.argv[1]])

	print("Analisando resultados nhmmer...")

	# Using readlines()
	file  = open('tmp.txt', 'r')
	lines = file.readlines()

	#Obteniendo listados de posibles mlts
	mlst = bestMLST(lines)

	output(mlst, sys.argv[1])


	

	print(f'acnb [{len(mlst[0])}]:\t{mlst[0][0]}')
	print(f'gap [{len(mlst[1])}]:\t{mlst[1][0]}')
	print(f'glta [{len(mlst[2])}]:\t{mlst[2][0]}')
	print(f'gyrb [{len(mlst[3])}]:\t{mlst[3][0]}')
	print(f'pgi [{len(mlst[4])}]:\t{mlst[4][0]}')
	print(f'rpoD [{len(mlst[5])}]:\t{mlst[5][0]}')

	os.remove("tmp.txt") 





def output( mlst, fastaFile):

	records = []
	basename = os.path.basename(fastaFile).split('.')

	for gen in mlst:
		for seq_record in SeqIO.parse( fastaFile, 'fasta'):
			if(seq_record.id == gen[0]):
				desc = seq_record.description.split(seq_record.id)
				rec  = SeqRecord( Seq(str(seq_record.seq)), id=seq_record.id, description= desc[1] )
				records.append(rec)

	SeqIO.write(records, f'{basename[0]}_MLST.fna', "fasta")




def bestMLST( lines ):
	acnb = []
	gap  = []
	glta = []
	gyrb = []
	pgi  = []
	rpod = []
	data = []
	
	for line in lines:
		if line.strip()[0] == '#':
			continue

		x = line.split()

		if x[2] == 'acnb':
			acnb.append(x[0])

		elif x[2] == 'gap':
			gap.append(x[0])

		elif x[2] == 'glta':
			glta.append(x[0])

		elif x[2] == 'gyrB':
			gyrb.append(x[0])

		elif x[2] == 'pgi':
			pgi.append(x[0])

		elif x[2] == 'rpoD':
			rpod.append(x[0])

		else:
			print("No se encontro MLTS")

	data.append(acnb)
	data.append(gap)
	data.append(glta)
	data.append(gyrb)
	data.append(pgi)
	data.append(rpod)

	return data


if __name__ == '__main__':
	main()