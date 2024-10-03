"""Script to produce shape frequency plot
The script takes three arguments, rna length, rna number to sample, and path to save plot
"""
import sys
# Needed for ssh to server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))

from src.foldingAlg import No_beam_list, BeamSearchDefault, Default_list
from src.framework import *

# torun = No_beam_list
torun = Default_list

if __name__ == "__main__":
    args = sys.argv
    print(*[x.name for x in torun])
    analysis = StrDensityAnalysis(100, 10, *torun)
    for i in range(int(args[1])):
        print("Running SDS", i, end='\r')
        analysis()
    df = analysis.to_df()
    df.to_csv(args[2]+'.csv', sep='\t', index=False)
    # analysis.plot_2d()
    # plt.savefig(args[2]+'.pdf', dpi=200, bbox_inches='tight')

# python script/shapefreq_wo_beam.py <nSamples> <outputfilenamebase>


