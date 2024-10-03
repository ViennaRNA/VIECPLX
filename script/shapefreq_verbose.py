"""Script to produce shape frequency plot
The script takes three arguments, rna length, rna number to sample, and path to save plot
"""
import sys
# Needed for ssh to server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from multiprocessing import Pool

sys.path.append(str(Path(__file__).parent.parent))

from src.foldingAlg import Default_list, BasicCoFold, BasicCoFoldRNAFold, LookBehindFold, BestHelixFoldRNAFold, BestHelixCoFold, FoldingRule, BeamSearchDefault

import RNA


args = sys.argv
seq_length = int(args[1])

OUT = Path(args[3])
OUT.mkdir(exist_ok=True)
running = [BeamSearchDefault]
# running = [BasicCoFold, BasicCoFoldRNAFold, BestHelixCoFold, BestHelixFoldRNAFold]
# RESULTS = [str(Path(OUT / "folding_result_{}.csv".format(x.name))) for x in Default_list]
RESULTS = [str(Path(OUT / "folding_result_{}.csv".format(x.name))) for x in running]

def read_sequences():
    with Path(OUT / "folding_result_RNAfold.csv").open() as f:
        for line in f.readlines():
            yield line.split('\t')[0]

# class Call:
#     def __init__(self):
#         self.foldings = [x for x in running]
#         self.results = [x for x in RESULTS]

#     def __call__(self, seq):
#         for folding, output in zip(self.foldings, self.results):
#             with open(output, 'a') as f:
#                 print(seq)
#                 print(seq, RNA.db_to_tree_string(folding.fold(seq), RNA.STRUCTURE_TREE_SHAPIRO_SHORT), sep='\t', file=f)
#         return None


def run(seq, ind):
    # print(ind)
    # seq = RNA.random_string(seq_length, "ACGU")
    print(seq, ind, end='\r')
    for folding, output in zip(running, RESULTS):
        with open(output, 'a') as f:
            dbn = folding.fold(seq)
            print(seq, RNA.db_to_tree_string(dbn, RNA.STRUCTURE_TREE_SHAPIRO_SHORT), dbn, sep='\t', file=f)
    return None


# with Pool(processes=8) as pool:
#     # pool.imap_unordered(run, range(int(args[2])))
#     pool.imap(Call(), (RNA.random_string(seq_length, 'ACGU') for _ in range(int(args[2]))))

# for ind in range(int(args[2])):
ind = 0
for seq in read_sequences():
    # seq = RNA.random_string(seq_length, 'ACGU')
    run(seq, ind)
    ind += 1
