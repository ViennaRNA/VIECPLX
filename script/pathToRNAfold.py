"""Script to produce shape frequency plot
The script takes three arguments, rna length, rna number to sample, and path to save plot
"""
import sys
from pathlib import Path
import dataclasses

import RNA

sys.path.append(str(Path(__file__).parent.parent))

# from src.foldingAlg import No_beam_list, Default_list
from src.foldingAlg import BeamSearchDefault, BestHelixFoldRNAFold, Final_list
from src.framework.neutral import PathToRNAfold

torun = Final_list[-1:0:-1]
# torun = No_beam_list

if __name__ == "__main__":
    args = sys.argv
    length = int(args[1])
    nb_seq = int(args[2])
    nb_rounds = int(args[3])
    sequences = [RNA.random_string(length, "ACGU") for _ in range(nb_seq)]
    # ExpObj = {folding.name: PathToRNAfold(folding, flag=RNA.STRUCTURE_TREE_SHAPIRO_SHORT) for folding in torun}
    # PathObj= PathToRNAfold(BestHelixFoldRNAFold)
    for folding in torun:
        # PathObj = PathToRNAfold(folding, flag=RNA.STRUCTURE_TREE_SHAPIRO_SHORT)
        PathObj = PathToRNAfold(folding)
        for seq in sequences:
            res = PathObj(seq, rounds=nb_rounds, store=False, loose=True)
            print(folding.name, *dataclasses.astuple(res), sep='\t', flush=True)
    # args = sys.argv1
    # analysis = Analysis(*torun)
    # LOG = f'{args[1]}_log.csv'
    # for ind in range(10):
    #     ROUND_PATH = f'{args[1]}_round{ind+1}.csv'
    #     upper = analysis.run_path_upper(100, 1, 50, nb_rounds=20, log_path=LOG)
    #     df = upper.to_df()
    #     df.to_csv(ROUND_PATH, sep='\t', index=False)

