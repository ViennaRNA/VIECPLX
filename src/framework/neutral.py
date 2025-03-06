"""Implementation of neutral network analysis
Schuster, P., Fontana, W., Stadler, P. F., & Hofacker, I. L. (1994). From sequences to shapes and back: a case study in RNA secondary structures. Proceedings of the Royal Society of London. Series B: Biological Sciences, 255(1344), 279-284.
"""

from collections import Counter
import random
import  dataclasses

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd
import seaborn as sns

import RNA

from .helper import cc_from_dbn, mutateSeq, Mutation, hamming_dist

MAXSTAY = 20


@dataclasses.dataclass
class FinalState:
    """Class object of final state of random walk
    """
    trial_seq: str
    seq: str
    dist: float
    total_moves: int
    total_tries: int
    hamming: int = dataclasses.field(init=False)

    def __post_init__(self):
        self.hamming = hamming_dist(self.trial_seq, self.seq)


class ShapeFreq:
    """Class object for shape frequency analysis
    """
    def __init__(self, folding, flag=RNA.STRUCTURE_TREE_SHAPIRO_SHORT):
        """Initiation with given folding algorithm

        Args:
            floding: flding algorithm
            flag: flag for coarse-grained shape
        """
        self.folding = folding
        self.shapes = Counter()
        self.flag = flag

    def __call__(self, seq):
        """Update shape counter for given sequence
        Given a sequence, a shape is defined as its Shapiro coarse grained folding, i.e. RNA.STRUCTURE_TREE_SHAPIRO_SHORT

        Args:
            seq: RNA sequence of interest
        """
        shape = RNA.db_to_tree_string(self.folding.fold(seq), self.flag)
        self.shapes.update([shape])
        return shape

    def rank(self):
        """Return shape ranks in desc order
        """
        return [x[1] for x in self.shapes.most_common()]


    def freq(self):
        """Return shape frequency in desc order
        """
        t = self.shapes.total()
        return [x[1]/t for x in self.shapes.most_common()]


    def rank_plot(self):
        """Plot each shape frequency in the order of frequency

        Sth like fig2 in paper
        """
        t = self.shapes.total()
        freq = [x[1]/t for x in self.shapes.most_common()]
        plt.plot(range(1, len(freq)+1), freq)


    def to_df(self):
        t = self.shapes.total()
        lst = list(zip(*[(x[0], x[1]/t) for x in self.shapes.most_common()]))
        return pd.DataFrame({'shape': lst[0], 'frequency': lst[1]})


def tree_edit_dist(s1, s2, flag=RNA.STRUCTURE_TREE_EXPANDED):
    """Return tree edit distance for given two structures in dbn
    This function converts first dbn to tree with flag, then call RNA.tree_edit_distance

    Args:
        s1: first structure in dbn
        s2: second structure in dbn
    """
    tree1 = RNA.make_tree(RNA.db_to_tree_string(s1, flag))
    tree2 = RNA.make_tree(RNA.db_to_tree_string(s2, flag))
    return int(RNA.tree_edit_distance(tree1, tree2))


class StrDensity:
    """Class object for structure density surface plot for single folding algorithm

    Args:
        folding: folding function
        ref: reference sequence
        str_distance: structure distance function
    """
    def __init__(self, folding, length):
        self.folding = folding
        self.ref_seq = None
        self.ref_str = None
        self.length = length
        self.bpdist= np.zeros((length+1, 2*length+1))
        self.treedist = np.zeros((length+1, 2*length+1))


    def __call__(self, seq, hamming=None):
        """Update distance result for given sequence

        Compute hamming distance to reference sequence if not given
        """
        if self.ref_seq is None:
            raise Exception('Reference sequence is not set yet')
        if hamming is None:
            hamming = sum(x != y for x, y in zip(self.ref_seq, seq))
        else:
            assert hamming == sum(x != y for x, y in zip(self.ref_seq, seq))

        dbn = self.folding.fold(seq)
        self.treedist[hamming, tree_edit_dist(self.ref_str, dbn)] += 1
        self.bpdist[hamming, RNA.bp_distance(self.ref_str, dbn)] += 1


    def set_ref(self, ref):
        """Set reference seqeunce
        """
        if not len(ref) == self.length:
            raise ValueError('Reference sequence length is different then initial one ({} vs {})'.format(len(ref), self.length))

        self.ref_seq = ref
        self.ref_str = self.folding.fold(ref)


    def to_df(self):
        """Export distances result as pandas.DataFrame
        """
        tmp = {'seq': [], 'str': [], 'bpcount': [], 'treecount': []}
        for i in range(self.length+1):
            for j in range(2*self.length+1):
                tmp['seq'].append(i)
                tmp['str'].append(j)
                tmp['bpcount'].append(self.bpdist[i,j])
                tmp['treecount'].append(self.treedist[i,j])
        return pd.DataFrame(tmp)


    def max_str_dist(self, label='bp'):
        """Return maximum structure distance
        """
        if label == 'bp':
            dist = self.bpdist
        else:
            dist = self.treedist
        return max(i for i in range(dist.shape[1]) if any(dist[:,i]))


    def plot_2d(self, label='bp'):
        if label == 'bp':
            count = 'bpcount'
        else:
            count = 'treecount'
        df = self.to_df()
        g = sns.jointplot(data=df, x='str', y='seq', kind='kde', xlim=(-1, self.max_str_dist()//10*10+1), ylim=(-1, 101), joint_kws={'weights': count, 'data': df})
        g.set_axis_labels('Structure distance', 'Sequence distance')
        g.ax_marg_x.remove()
        g.ax_marg_y.remove()
        return g


    def plot_3d(self, figsize=(10,10), label='bp', **kwargs):
        """Draw SDS using pyplot 3d surface
        keywords are given to plot_surface
        """
        # x for sequence distance
        # y for structure distance
        if label == 'bp':
            dist = self.bpdist
        else:
            dist = self.treedist
        ymin = 0
        ymax = self.max_str_dist() // 10 * 10 + 1
        xmin = 0
        xmax = self.length+1
        xx, yy = np.mgrid[xmin:xmax, ymin:ymax]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        total = np.sum(dist)
        f = np.reshape([dist[i, j]/total for i, j in zip(*positions)], xx.shape)
        f[np.where(f==0)] = np.nan

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=figsize)
        surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, alpha=.5)
        ax.set_ylabel('Structure distance')
        ax.set_xlabel('Sequence distance')
        ax.set_zlabel('Density')
        ax.view_init(15, -25)
        ax.set_xlim(ax.get_xlim()[::-1])
        return fig

class AdaptiveWalk:
    """Class object for general adaptive walk in sequence space
    """
    def __init__(self, folding, debug=False, flag=None):
        self.folding = folding
        self.result = []
        self.debug = debug
        self.flag = flag
        self.cur_state = {}
        self._mutation_bases = None

    def _max_dist(self, trial_seq):
        """Maximum distance given trial_seq, default to length
        The value doesn't need to precise
        """
        return len(trial_seq)

    def convert_str(self, dbn):
        """Convert structure according to flag
        """
        if self.flag is not None:
            return RNA.db_to_tree_string(dbn, self.flag)
        return dbn

    def _is_better(self, best_state, new_state):
        """Return True if new_state is better
        """
        return new_state.dist < best_state.dist

    def _stop_rounds(self, best_state):
        return best_state.dist == 0

    def __call__(self, trial_seq, rounds=1, loose=False, store=True):
        """Compute path for given target and trial sequence with several rounds
        """
        # Store the best state among rounds
        best = None
        for _ in range(rounds):
            to_call = self.singleCallLoose if loose else self.singleCall
            final = to_call(trial_seq)
            if best is None or self._is_better(best, final):
                best = final
            # Stop the round if the reached state is the possible best state
            if self._stop_rounds(best):
                break

        if store:
            self.result.append(best)
        return best

    def distance(self, cur_seq):
        """Compute distance of current sequence
        """
        return 0

    def all_possible_moves(self):
        """Create next moves based on current state
        """
        cur_seq = self.cur_state['seq']
        if self._mutation_bases is None:
            self._mutation_bases = [(x, y) for x in range(len(cur_seq)) for y in range(3)]
        random.shuffle(self._mutation_bases)
        for pos, offset in self._mutation_bases:
            new_seq = cur_seq[:pos] + Mutation[cur_seq[pos]][offset] + cur_seq[pos+1:]
            yield new_seq


    def nb_possible_moves(self, trial_seq):
        # return 3 * len(trial_seq)
        return len(trial_seq)


    def is_valid(self, seq):
        return True

    def init_state(self, trial_seq):
        """Create init state
        """
        # Clean up current state
        self.cur_state = {}
        self.cur_state['total_moves'] = 0
        self.cur_state['total_tries'] = 0
        self.cur_state['trial_seq'] = trial_seq
        self.cur_state['seq'] = trial_seq
        self.cur_state['dist'] = self.distance(trial_seq)

    def singleCall(self, trial_seq):
        """Compute path for given target and trial sequence with one trial only
        """
        self.init_state(trial_seq)

        total_moves = 0
        total_tries = 0
        stop = self.cur_state['dist'] == 0

        while (not stop):
            stop = True
            for new_seq in self.all_possible_moves():
                total_tries += 1
                new_dist = self.distance(new_seq)
                # Check whether the new move is valid
                if new_dist >= self.cur_state['dist'] or not self.is_valid(new_seq):
                    continue
                # Accept move when distance strickly decreases
                # ignore the rest of moves
                if new_dist < self.cur_state['dist']:
                    self.cur_state['dist'] = new_dist
                    self.cur_state['seq'] = new_seq
                    total_moves += 1
                    # Keep walking as distance is not yet null
                    self.cur_state['total_moves'] = total_moves
                    self.cur_state['total_tries'] = total_tries
                    if self.cur_state['dist'] > 0:
                        stop = False
                    break
        return FinalState(**{field.name: self.cur_state[field.name] for field in dataclasses.fields(FinalState) if field.init})


    def singleCallLoose(self, trial_seq, max_stay=20):
        """Loose version of adaptive walk. The neutral move is allowed, i.e. move that does not improve distance.
        The neutral move is forbidden when enough neutral moves (max_stay) are made
        """
        self.init_state(trial_seq)

        nb_possible_moves = self.nb_possible_moves(trial_seq)
        total_moves = 0
        total_tries = 0
        stop = self.cur_state['dist'] == 0
        stay = 0


        while (not stop):
            stop = True
            for ind, new_seq in enumerate(self.all_possible_moves()):
                total_tries += 1
                # Check whether the new move is valid
                new_dist = self.distance(new_seq)
                cur_dist = self.cur_state['dist']
                if new_dist > cur_dist or (new_dist == cur_dist and stay >= MAXSTAY) or (not self.is_valid(new_seq)):
                    continue

                # Neutral move, and cummulated stay number is smaller than MAXSTAY
                # Conditionally accept neutral move (distance stays the same)
                elif new_dist == cur_dist:
                    # We don't want to accept neutral move at early stage
                    if random.uniform(0, 1) < (ind+1)/nb_possible_moves:
                        stay += 1
                    else:
                        continue

                # The move improves distance
                else:
                    stay = 0
                # Here, we accept the new move
                self.cur_state['dist'] = new_dist
                self.cur_state['seq'] = new_seq
                total_moves += 1
                self.cur_state['total_moves'] = total_moves
                self.cur_state['total_tries'] = total_tries
                # Keep walking as distance is not yet null
                if self.cur_state['dist'] > 0:
                    stop = False
                break

        return FinalState(**{field.name: self.cur_state[field.name] for field in dataclasses.fields(FinalState) if field.init})

class AdaptiveWalkStructure(AdaptiveWalk):
    """Class object of adaptive walk with a structure as reference
    """
    def set_reference(self, structure):
        self._ref_ss = structure
        self._cc_of_ref = cc_from_dbn(self._ref_ss)

    def __call__(self, trial_seq, ref_ss, rounds=1, loose=False, store=True):
        self.set_reference(ref_ss)
        return super().__call__(trial_seq, rounds=rounds, loose=loose, store=store)

    def is_valid(self, seq):
        """Make sure sequence is in the neutral network
        """
        return self.folding.fold(seq) == self._ref_ss

    def all_possible_moves(self):
        """Generate all possible sequence w.r.t the reference structure
        """
        random.shuffle(self._cc_of_ref)
        for to_mutate in self._cc_of_ref:
            if len(to_mutate) == 1:
                possible = list('ACGU')
            elif len(to_mutate) == 2:
                possible = ['AU', 'CG', 'GC', 'GU', 'UA', 'UG']
            else:
                raise Exception("Mutation positions is more than two")
            random.shuffle(possible)
            for x in possible:
                if x == ''.join(self.cur_state['seq'][i] for i in to_mutate):
                    continue
                new_seq = self.cur_state['seq']
                for i, c in zip(to_mutate, x):
                    new_seq = new_seq[:i] + c + new_seq[i+1:]
                yield new_seq


class PathUpperLegacy:
    """Class object to compute upper bound for neutral path
    """
    def __init__(self, folding, debug=False, flag=None):
        self.folding = folding
        self.path = []
        self.debug = debug
        self.flag = flag


    def convert_str(self, dbn):
        """Convert structure according to flag
        """
        if self.flag is not None:
            return RNA.db_to_tree_string(dbn, self.flag)
        return dbn

    def __call__(self, target_seq, trial_seq, rounds=1, loose=False):
        """Compute path for given target and trial sequence with several rounds
        """
        current = len(target_seq)
        for _ in range(rounds):
            if loose:
                upper = self.singleCallLoose(target_seq, trial_seq)
            else:
                upper = self.singleCall(target_seq, trial_seq)
            current = min(current, upper)
        self.path.append(current)


    def singleCall(self, target_seq, trial_seq):
        """Compute path for given target and trial sequence with one trial only
        """
        ref_str = self.folding.fold(trial_seq)
        cc_ref = cc_from_dbn(ref_str)
        if self.debug:
            print(target_seq)
            print(ref_str)
            print(trial_seq)

        # Here we create list to mutate
        mutate_list = []
        for t in cc_ref:
            if not all(target_seq[x] == trial_seq[x] for x in t):
                mutate_list.append([(x, target_seq[x]) for x in t])

        cur_trial_seq = trial_seq
        stop = False

        while (not stop):
            stop = True
            for ind in np.random.permutation(len(mutate_list)):
                # New_sequence is closer to target sequence
                new_seq = cur_trial_seq
                for x, c in mutate_list[ind]:
                    new_seq = new_seq[:x] + c + new_seq[x+1:]
                # Set new trial and keep running
                if self.convert_str(ref_str) == self.convert_str(self.folding.fold(new_seq)):
                    if self.debug:
                        print(new_seq)
                    del mutate_list[ind]
                    cur_trial_seq = new_seq
                    stop = False
                    break
        return sum(x!=y for x, y in zip(target_seq, cur_trial_seq))


    def singleCallLoose(self, target_seq, trial_seq):
        """Compute path for given target and trial sequence with one trial only while allowing same hamming distance new trial in the walk
        This version follows the same implementation as original one
        """
        ref_str = self.folding.fold(trial_seq)
        cc_ref = cc_from_dbn(ref_str)
        n_pos = len(cc_ref)
        if self.debug:
            print(target_seq)
            print(ref_str)
            print(trial_seq)


        cur_trial_seq = trial_seq
        stop = False
        stay = 0


        while (not stop):
            stop = True
            dist = None
            p_dist = hamming_dist(target_seq, cur_trial_seq)
            for l, ind in enumerate(np.random.permutation(n_pos)):
                # Position(s) to mutate
                lst = cc_ref[ind]
                trial = ''.join(cur_trial_seq[x] for x in lst)
                target = ''.join(target_seq[x] for x in lst)
                if trial == target:
                    continue
                # New_sequence is closer to target sequence
                for mutations, d in possible_mutations(trial, target):
                    # Mutate current trial
                    new_seq = cur_trial_seq
                    for x, c in zip(lst, mutations):
                        new_seq = new_seq[:x] + c + new_seq[x+1:]
                    # Conditionally accept no improve mutation
                    if d == 0:
                        if stay < MAXSTAY and random.uniform(0, 1) < (l+1)/n_pos:
                            stay += 1
                        else:
                            continue

                    # Set new trial and keep running if structure of new_trial is same as the ref
                    if self.convert_str(ref_str) == self.convert_str(self.folding.fold(new_seq)):
                        if self.debug:
                            print(new_seq)
                        cur_trial_seq = new_seq
                        stop = False
                        dist = d
                        break

                if not stop:
                    n_dist = hamming(target_seq, cur_trial_seq)
                    assert n_dist - p_dist == dist
                    if dist < 0:
                        stay = 0

                    break
        return hamming_dist(target_seq, cur_trial_seq)

class PathUpper(AdaptiveWalkStructure):
    """Class object to compute upper bound for neutral path
    """
    def distance(self, cur_seq):
        """Distance to the target sequence
        """
        return hamming_dist(cur_seq, self._target_seq)

    def __call__(self, target_seq, trial_seq, rounds=1, loose=False, store=True):
        """Compute path for given target and trial sequence with several rounds
        """
        self._target_seq = target_seq
        ref_ss = self.folding.fold(trial_seq)
        return super().__call__(trial_seq, ref_ss, rounds=rounds, loose=loose, store=store)

class PathUpperRNAfold(PathUpper):
    """Class object to compute upper bound for neutral path
    """
    def is_valid(self, seq):
        """Make sure sequence is in the neutral network and same prediction as RNAfold
        """
        # First check whether is in the same neutral network
        # Next check whether is also RNAfold prediction
        flag = (self.folding.fold(seq) == self._ref_ss)
        if flag:
            fc = RNA.fold_compound(seq)
            mfe_ss, mfe = fc.mfe()
            return mfe == fc.eval_structure(self._ref_ss)
        return False


class PathToRNAfold(AdaptiveWalk):
    """Class object for a adaptive walk within sequence space till given folding algorithm predicts same as RNAfold
    """
    def init_state(self, trial_seq):
        super().init_state(trial_seq)

    def distance(self, cur_seq):
        if self.flag is None:
            return self.distance_bp(cur_seq)
        return self.distance_tree(cur_seq)

    def distance_tree(self, cur_seq):
        fc = RNA.fold_compound(cur_seq)
        ss = self.folding.fold(cur_seq)
        tree_ss = RNA.make_tree(RNA.db_to_tree_string(ss, self.flag))

        dist = []
        for x in fc.subopt(0):
            tree_mfe = RNA.make_tree(RNA.db_to_tree_string(x.structure, self.flag))
            dist.append(RNA.tree_edit_distance(tree_mfe, tree_ss))
        return min(dist)

    def distance_bp(self, cur_seq):
        """BP distance between stuructures of cur_seq computed with folding algorithm and RNAfold
        We further check whether the one predicted by folding algorithm is one of the MFE
        """
        fc = RNA.fold_compound(cur_seq)
        mfe_ss, mfe = fc.mfe()
        ss = self.folding.fold(cur_seq)
        dist = RNA.bp_distance(ss, mfe_ss)
        # If the structure is distinct with RNAfold, check whether the structure is one of the MFE
        if dist > 0 and mfe < fc.eval_structure(ss):
            return dist
        return 0

    def _is_better(self, best_state, new_state):
        """New state is better if distance is smaller or hamming distance is smaller in case of equal distance
        """
        return (new_state.dist < best_state.dist) or (new_state.dist == best_state.dist and new_state.hamming < best_state.hamming)

    def _stop_rounds(self, best_state):
        return best_state.dist == 0 and best_state.hamming <= 1


class MaxNeutralPath(AdaptiveWalkStructure):
    """Class object for starting from a sequence goes as far as possible within neutral network
    """
    def distance(self, cur_seq):
        """Distance is defined as the amount of remaining mutations needed to reach full change from trial sequence
        """
        return len(cur_seq) - hamming_dist(cur_seq, self.cur_state['trial_seq'])

    def __call__(self, trial_seq, rounds=1, loose=False, store=True):
        ref_ss = self.folding.fold(trial_seq)
        return super().__call__(trial_seq, ref_ss, rounds=rounds, loose=loose, store=store)


class MaxNeutralPathWithRNAfold(MaxNeutralPath):
    """Class object for starting from a sequence goes as far as possible within neutral network
    We ensure during the all path, the target structure is RNAfold MFE
    """
    def is_valid(self, seq):
        """Make sure sequence is in the neutral network and same prediction as RNAfold
        """
        # First check whether is in the same neutral network
        # Next check whether is also RNAfold prediction
        flag = (self.folding.fold(seq) == self._ref_ss)
        if flag:
            fc = RNA.fold_compound(seq)
            mfe_ss, mfe = fc.mfe()
            return mfe == fc.eval_structure(self._ref_ss)
        return False


    def __call__(self, trial_seq, rounds=1, loose=False, store=True):
        # Here we first mutate trial seq such that the prediction is same as RNAfold
        tmp = PathToRNAfold(self.folding)
        final = tmp(trial_seq, rounds=10, loose=True)
        # Cannot reach to the same prediction
        if final.dist > 0:
            return None
        return super().__call__(final.seq, rounds=rounds, loose=loose, store=store)


class NeighborHood:
    """Class object to compute the neighborhood w.r.t. RNAfold
    """
    def __init__(self, folding, debug=False, flag=None):
        self.folding = folding
        self.path = []
        self.debug = debug
        self.flag = flag
        self.dist = None


    def convert_str(self, dbn):
        """Convert structure according to flag
        """
        if self.flag is not None:
            return RNA.db_to_tree_string(dbn, self.flag)
        return dbn

    def __call__(self, target_seq, trial_seq, rounds=1, loose=False):
        """Compute path for given target and trial sequence with several rounds
        """
        current = len(target_seq)
        for _ in range(rounds):
            if loose:
                upper = self.singleCallLoose(target_seq, trial_seq)
            else:
                upper = self.singleCall(target_seq, trial_seq)
            current = min(current, upper)
        self.path.append(current)

    # TODO: implement shapre comparaison
    def assert_equal(ss1, ss2):
        """Compare two structure
        """
        return ss1 == ss2

    def singleCall(self, trial_seq, sample_size=10000):
        """Design the target structure from given start_seq
        """
        # Assert folding algorithm has the same prediction as RNAfold on trial_seq
        ss1 = self.folding.fold(trial_seq)
        ss2 = RNA.fold(trial_seq)

        if not assert_equal(ss1, ss2):
            raise ValueError(f"The predicted structure for the input sequence is different than the one of RNAfold\n{ss1} {self.folding.name}\n{ss2} RNAfold")

        same_count = 0

        for _ in range(sample_size):
            seq = mutateSeq(trial_seq)
            same_count += self.assert_equal(self.folding.fold(seq), RNA.fold(seq))
        return same_count

def possible_mutations(trial, target):
    """Generate possible mutation list
    """
    if len(trial) == 1:
        possible = list('ACGU')
    elif len(trial) == 2:
        possible = ['AU', 'CG', 'GC', 'GU', 'UA', 'UG']
    else:
        raise Exception("Mutation positions is more than two")
    random.shuffle(possible)
    for x in possible:
        # Mutation should be different than the origin one
        if x == trial:
            continue
        d = distChange(target, trial, x)
        # Only allow mutation w/o increasing distance
        if d <= 0:
            yield x, d


def distChange(target, old_trial, new_trial):
    """Return the change of distance to the target from old to new trial
    """
    return hamming_dist(new_trial, target) - hamming_dist(old_trial, target)

