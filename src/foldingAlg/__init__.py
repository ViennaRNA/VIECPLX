
from . import basic
from . import basic_co_folding
from . import best_helix_co_folding
from . import RNA_folding_rule
from . import look_behind
from . import beam_search_fold


from .basic import *


# Here we define all folding algorithms

class BeamSearch(Folding):
    def __init__(self, name, bp_threshold=0.01, l=50, bps=-1):
        super().__init__(name, None)
        self.beam = beam_search_fold.BeamSearch()
        self.bp_threshold=bp_threshold
        self.l = l
        self.bps = bps
 
    def fold(self, w):
        return self.beam.fold(w, bp_threshold=self.bp_threshold, l=self.l, bps=self.bps, debug=False)


class CoFoldRNA(Folding):
    def fold(self, w):
        if self._fold is not None:
            return self._fold(w, helix_length=6, useRNAfoldForRound2=True)
        else:
            raise TypeError('Folding algorithm is undefined')


BasicCoFold          = Folding('Basic cofold', basic_co_folding.fold)
BasicCoFoldRNAFold   = CoFoldRNA('Basic cofold with RNAfold in second step', basic_co_folding.fold)
BestHelixCoFold      = Folding('Best helix cofold', best_helix_co_folding.fold)
BestHelixFoldRNAFold = CoFoldRNA('Best helix cofold with RNAfold in second step', best_helix_co_folding.fold)
FoldingRule          = Folding('Folding rule', RNA_folding_rule.fold)
LookBehindFold       = Folding('Look behind fold', look_behind.fold)
BeamSearchDefault    = BeamSearch('Beam search')
BeamSearchLegacy     = BeamSearch('Beam search legacy')



Default_list = [ViennaFold, BasicCoFold, BestHelixCoFold, FoldingRule, LookBehindFold, BasicCoFoldRNAFold, BestHelixFoldRNAFold, BeamSearchDefault]
No_beam_list = [ViennaFold, BasicCoFold, BestHelixCoFold, FoldingRule, LookBehindFold, BasicCoFoldRNAFold, BestHelixFoldRNAFold]
Final_list = [ViennaFold, LookBehindFold, BasicCoFold, BestHelixFoldRNAFold, FoldingRule, BeamSearchDefault]
