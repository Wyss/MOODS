import unittest
import sys
from os.path import join, abspath, dirname

LOCAL_DIR = abspath(dirname(__file__))
sys.path.append(abspath(join(LOCAL_DIR, '..')))

import MOODS

import fasta

DIST_DIR = abspath(dirname(dirname(LOCAL_DIR)))
print(DIST_DIR)
fasta_filepath = join(DIST_DIR, "examples/data/sequence/dnaACGT.txt")
records = fasta.parseFasta(fasta_filepath)

seq = records[0][1]

matrix1 = [     [0,1,0,0,0,0,0,1,1,0],
                [1,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0],
                [0,0,1,1,1,1,1,0,0,1]
            ]
matrix2 = [     [10,0,10,3,5,5],
                [0,5,0,3,5,0,5],
                [0,1,0,3,0,5,0],
                [0,4,0,1,0,0,5]
            ]

results = MOODS.search(seq, [matrix1, matrix2], 0.011)

print("Matrix 1 results: "+ str(len(results[0])))
print("Matrix 2 results: "+ str(len(results[1])))


matrices = [matrix1, matrix2]
thresholds = [0.011, 0.011]
bg = MOODS.bg_from_sequence(seq, 0.1)
q = 7
absolute_threshold = False
both_strands=False
ms = MOODS.MOODSSearch(matrices, thresholds, bg, q, absolute_threshold, both_strands)
results = ms.search(seq)

print("New Matrix 1 results: "+ str(len(results[0])))
print("New Matrix 2 results: "+ str(len(results[1])))

ms2 = MOODS.MOODSSearch(matrices, thresholds, bg, q, absolute_threshold, both_strands)