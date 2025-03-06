"""Script to produce shape frequency plot
The script takes three arguments, rna length, rna number to sample, and path to save plot
"""
import sys
from pathlib import Path
import dataclasses

import RNA

sys.path.append(str(Path(__file__).parent.parent))

# from src.foldingAlg import No_beam_list, Default_list
from src.foldingAlg import BeamSearchDefault, BestHelixFoldRNAFold, ViennaFold, FoldingRule, LookBehindFold, BasicCoFold
from src.framework.neutral import MaxNeutralPathWithRNAfold

# torun = Final_list[-1::-1]
# torun = Final_list
torun = LookBehindFold
# torun = No_beam_list

if __name__ == "__main__":
    args = sys.argv
    length = int(args[1])
    nb_seq = int(args[2])
    nb_rounds = int(args[3])
    # sequences = [RNA.random_string(length, "ACGU") for _ in range(nb_seq)]
    PathObj = MaxNeutralPathWithRNAfold(torun)
    ind = 0
    while ind < nb_seq:
        seq = RNA.random_string(length, "ACGU")
        res = PathObj(seq, rounds=nb_rounds, store=False, loose=True)
        if res is None:
            continue
        print(torun.name, *dataclasses.astuple(res), sep='\t', flush=True)
        ind +=1
