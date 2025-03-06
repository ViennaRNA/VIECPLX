"""Script to produce shape frequency plot
The script takes three arguments, rna length, rna number to sample, and path to save plot
"""
import sys
# Needed for ssh to server
from pathlib import Path
import dataclasses

import RNA

sys.path.append(str(Path(__file__).parent.parent))

from src.foldingAlg import BeamSearchDefault, BestHelixFoldRNAFold, ViennaFold, FoldingRule, LookBehindFold, BasicCoFold
from src.framework.neutral import PathUpperRNAfold, PathToRNAfold

# torun = Default_list
torun = BeamSearchDefault
# torun = [ViennaFold]
# torun = No_beam_list

if __name__ == "__main__":
    args = sys.argv
    length = int(args[1])
    nb_seq = int(args[2])
    nb_rounds = int(args[3])
    # sequences = [RNA.random_string(length, "ACGU") for _ in range(nb_seq)]
    PathObj = PathUpperRNAfold(torun)
    tmp = PathToRNAfold(torun)
    while True:
        seq1 = RNA.random_string(length, "ACGU")
        final = tmp(seq1, rounds=10, loose=True)
        if final.dist > 0:
            continue
        target_seq = final.seq
        break

    ind = 0
    while ind < nb_seq:
        seq2 = RNA.random_string(length, "ACGU")
        # Cannot reach to the same prediction
        final = tmp(seq2, rounds=10, loose=True)
        # Cannot reach to the same prediction
        if final.dist > 0:
            continue
        trial_seq = final.seq
        res = PathObj(target_seq, trial_seq, rounds=nb_rounds, store=False, loose=True)
        print(torun.name, target_seq, *dataclasses.astuple(res), sep='\t', flush=True)
        ind +=1

