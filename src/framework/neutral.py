"""Implementation of neutral network analysis
Schuster, P., Fontana, W., Stadler, P. F., & Hofacker, I. L. (1994). From sequences to shapes and back: a case study in RNA secondary structures. Proceedings of the Royal Society of London. Series B: Biological Sciences, 255(1344), 279-284.
"""

from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd
import seaborn as sns

import RNA

from .helper import cc_from_dbn


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


class PathUpper:
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

    def __call__(self, target_seq, trial_seq, rounds=1):
        """Compute path for given target and trial sequence with several rounds
        """
        current = len(target_seq)
        for _ in range(rounds):
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

