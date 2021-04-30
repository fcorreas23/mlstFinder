#! /usr/bin/env python
import os, sys, argparse, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def getArgs():
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='<FILE>', help='Input sequence fasta file.')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-t', '--threads', type=int, metavar='int', default= 2, help='number of parallel CPU workers to use for multithreads  [2]')
    optional.add_argument('-e', '--evalue', type=float , metavar='float' , default = 10.0, help='report sequences <= this e-value threshold in output  [10.0]  (x>0)')
    optional.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<FILE>', help='Output file name. Default: Input file name + _hmm.txt')

    args = parser.parse_args()

    return args

def parseHmmerResult( tblout ):
    acnb = []
    gap  = []
    glta = []
    gyrb = []
    pgi  = []
    rpod = []
    data = []
    
    for line in tblout:
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
            print("No se encontro MLST")
    
    data.append(acnb)
    data.append(gap)
    data.append(glta)
    data.append(gyrb)
    data.append(pgi)
    data.append(rpod)
    
    return data


def makeFasta( mlst, fastaFile):
    genes = []
    sequences = []
    for gen in mlst: 
        for seq_record in SeqIO.parse( fastaFile, 'fasta'):
            if(seq_record.id == gen[0]):
                desc = seq_record.description.split(seq_record.id)
                rec  = SeqRecord( Seq(str(seq_record.seq)), id=seq_record.id, description= desc[1] )
                genes.append(rec)
                sequences.append(str(seq_record.seq))
    
    return genes, sequences


def main():
    #MAIN
    args = getArgs()
   
    filename = os.path.basename( args.input)
    basename = filename.split('.')
    print(f'\nBuscando genes MLST ( acnB, gap, glta, gyrB, pgi, rpoD ) en {filename}')
    print('Ejecutando nhmmer....')
    subprocess.call(['nhmmer','--cpu', str(args.threads), '-E', str(args.evalue) ,'-o', 'result_nhmmer.txt','--tblout', 'tmp.txt', 'hmm/mlts.hmm', args.input ])
    
    print("Analizando resultados nhmmer...")
    tblout = open('tmp.txt')
    lines = tblout.readlines()
    mlst = parseHmmerResult(lines)
    print("Generando archivos fasta MLST...")
    genes, sequences = makeFasta( mlst, args.input)
    SeqIO.write(genes, f'{basename[0]}_MLST.fna', "fasta")
    print("Concatenando archivos fasta MLST...")
    sequence = "".join(sequences)
    rec1 = SeqRecord( Seq(str(sequence)), id='S01_mlst', description= '[acnB, gap, glta, gyrB, pgi, rpoD]' )
    SeqIO.write(rec1, 'S01_MLST.fna', "fasta")
    os.remove("tmp.txt") 
    
    print('\nMLST[cant.]\tSeqID')
    print('-----------------------')
    print(f'acnb [{len(mlst[0])}]:\t{mlst[0][0]}')
    print(f'gap [{len(mlst[1])}]:\t{mlst[1][0]}')
    print(f'glta [{len(mlst[2])}]:\t{mlst[2][0]}')
    print(f'gyrb [{len(mlst[3])}]:\t{mlst[3][0]}')
    print(f'pgi [{len(mlst[4])}]:\t{mlst[4][0]}')
    print(f'rpoD [{len(mlst[5])}]:\t{mlst[5][0]}')

if __name__ == '__main__':
    main()