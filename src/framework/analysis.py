"""Script for general analysis
"""

from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import RNA

from .neutral import ShapeFreq, StrDensity, PathUpper
from .continuity import ContinuousEvol
from .helper import mutateSeq



# Here we define the function/feature that need to be equally applied to all folding algorithms
# For example, the initial sequence of pop evolution
class Analysis:
    def __init__(self, *folding):
        """Initiate analysis object for given folding objects
        """
        self.foldings = {f.name: f for f in folding}

    def run_shapes_frequency(self, length, nSamples, report_equal=False, **kwargs):
        """Compute shapes frequency in sequence space

        Args:
            length: sequence length of interest
            nSamples: number of sampled sequences
        """
        shapes = ShapeFreqAnalysis(length, *self.foldings.values(), **kwargs)
        for i in range(nSamples):
            print("Running Shape Freq", i, end='\r')
            res = shapes(report_equal=report_equal)
            if res is not None:
                print(res)
        return shapes


    def run_path_upper(self, length, nb_target, nb_trial):
        """Run neutral path upper bound analysis

        Args:
            length: length of RNA
            nb_target: number of target sequenes
            nb_trial: nuber of trial sequences for each target sequence
        """
        upper = PathUpperAnalysis(*self.foldings.values())
        for _ in range(nb_target):
            target_seq = RNA.random_string(length, 'ACGU')
            for _ in range(nb_trial):
                trial_seq = RNA.random_string(length, 'ACGU')
                upper(target_seq, trial_seq)
        return upper


class ShapeFreqAnalysis:
    def __init__(self, n, *folding, flag=RNA.STRUCTURE_TREE_SHAPIRO_SHORT, save_seq=False):
        self.records = {f.name: ShapeFreq(f, flag) for f in folding}
        self.length = n
        self.save_seq = save_seq
        self.sequences = []

    def __call__(self, report_equal=False):
        """Randomly sample one sequence then update result
        """
        seq = RNA.random_string(self.length, "ACGU")
        if self.save_seq:
            self.sequences.append(seq)
        shapes = [x(seq) for x in self.records.values()]
        fst = shapes[0]
        fst_I = fst.replace('B', 'I')
        if not report_equal:
            return None
        if all(x==fst for x in shapes[1:]):
            return seq, fst
        if all(x.replace('B', 'I')==fst_I for x in shapes[1:]):
            return seq, fst_I
        return None

    def to_df(self):
        """SDS result in pd.DataFrame
        """
        tmp = []
        for name, x in self.records.items():
            df = x.to_df()
            df.loc[:, 'label'] = name
            tmp.append(df)
        return pd.concat(tmp, ignore_index=True)


    def plot(self):
        for name, count in self.records.items():
            plt.plot(count.freq(), label=name)

        plt.legend()
        plt.xlabel('Rank')
        plt.ylabel('Frequency')
        plt.xscale('log')
        plt.yscale('log')


class StrDensityAnalysis:
    """Structure density surface analysis for several folding algorithms

    Args:
        nbClass: number of sampled sequences for each hamming class
        save_fig: flag to record all reference sequences
    """
    def __init__(self, length, nbClass, *folding, save_seq=False):
        self.records = {f.name: StrDensity(f, length) for f in folding}
        self.length = length
        self.nbClass = nbClass
        self.save_seq = save_seq
        self.ref = []


    def __call__(self):
        """Initiate one reference sequence and compute structure density surface
        """
        ref = RNA.random_string(self.length, "ACGU")
        if self.save_seq:
            self.ref.append(ref)
        for x in self.records.values():
            x.set_ref(ref)

        for hamming in range(1, self.length+1):
            for _ in range(self.nbClass):
                seq = mutateSeq(ref, hamming=hamming)
                for x in self.records.values():
                    x(seq, hamming=hamming)


    def to_df(self):
        """SDS result in pd.DataFrame
        """
        tmp = []
        for name, x in self.records.items():
            df = x.to_df()
            df.loc[:, 'label'] = name
            tmp.append(df)
        return pd.concat(tmp, ignore_index=True)



    def plot_2d(self):
        """Put SDS of all folding algorithms in one 2D plot
        """
        allDf = self.to_df()
        strlim = max(allDf['str'])
        g = sns.jointplot(data=allDf, x='str', y='seq', hue='label', kind='kde', xlim=(-1,strlim//10*10+1), ylim=(-1, 101), joint_kws={'weights': 'count', 'data': allDf})
        g.set_axis_labels('Structure distance', 'Sequence distance')
        g.ax_marg_x.remove()
        g.ax_marg_y.remove()
        g.ax_joint.get_legend().set_title(None)
        return g


class PathUpperAnalysis:
    def __init__(self, *folding):
        self.records = {f.name: PathUpper(f) for f in folding}
        self.target_seq = []
        self.trial_seq = []


    def __call__(self, target_seq, trial_seq):
        self.target_seq.append(target_seq)
        self.trial_seq.append(trial_seq)
        for x in self.records.values():
            x(target_seq, trial_seq)


    def to_df(self):
        tmp = []
        for name, x in self.records.items():
            tmp.append(pd.DataFrame({'label': name, 'path': x.path}))
        return pd.concat(tmp, ignore_index=True)

class FoldingHierarchy:
    """Hieracy of structure prediction algorithms

    Args:
        length: sequence length
        n: sequence sample size
        folding: folding algorithms to compare
    """
    def __init__(self, length, n, *folding):
        self.records = []
        self.length = length
        self.n = n
        self.folding = folding


    def sample(self, tree_rep = RNA.STRUCTURE_TREE_SHAPIRO_SHORT):
        self.records = []
        for _ in range(self.n):
            record = {}
            seq = RNA.random_string(self.length, "ACGU")
            record['seq'] = seq
            for f in self.folding:
                db = f.fold(seq)
                tree = RNA.db_to_tree_string(db,tree_rep)
                record[f'{f.name} db'] = db
                record[f'{f.name} tree'] = tree
                
                if  db == record['RNAfold db']: 
                    c_db = True
                else:
                    c_db = False
                record[f'{f.name} correctDB'] = c_db
                if  tree == record['RNAfold tree']: 
                    c_tree = True
                else:
                    c_tree = False
                record[f'{f.name} correctShape'] = c_tree
                
            self.records.append(record)
        return self.records

    def to_csv(self, out_path = 'hierarchy_records.tsv'):
        df = pd.DataFrame(self.records)
        df.to_csv(out_path , sep='\t')

    def plot_correct_vrna(self, figure_path = 'hierachy_correct'):
        df = pd.DataFrame(self.records)
        data = []
        for f in self.folding:
            data.append({'algorithm':f.name, 'correct db fraction': df[f'{f.name} correctDB'].sum()/self.n, 'correct shape fraction':df[f'{f.name} correctShape'].sum()/self.n})            
        df = pd.DataFrame(data)
        fig, ax = plt.subplots(1,1,figsize=(8, 5))
        df.plot.bar(x = 'algorithm', y = 'correct shape fraction', ax = ax)
        fig.savefig(figure_path+'.pdf', bbox_inches='tight')
        plt.yscale('log')
        fig.savefig(figure_path+'_log.pdf', bbox_inches='tight')
        return df

    def plot_probabilities(self, figure_path = 'hierachy_prob'):
        df = pd.DataFrame(self.records)
        data = []
        for f1 in self.folding:
            for f2 in self.folding:
                p_f1 = df[f'{f1.name} correctShape'].sum()/self.n
                p_f2 = df[f'{f2.name} correctShape'].sum()/self.n

                df_f2 = df.loc[df[f'{f2.name} correctShape']]
                p_f1_cond_f2 = df_f2[f'{f1.name} correctShape'].sum()/len(df_f2)

                data.append({'f1': f1.name, 'f2': f2.name, 'P(f1)': p_f1, 'P(f2)':p_f2, 'P(f1|f2)':p_f1_cond_f2 })
                
        df = pd.DataFrame(data)
        pivoted=pd.pivot_table(df, values="P(f1|f2)", index=['f1'], columns=["f2"])
        pivoted.to_csv(figure_path + '.tsv', sep='\t')

        fig, ax = plt.subplots(1,1,figsize=(5, 5))
        sns.heatmap(pivoted, square =True, cmap = 'Blues')
        fig.savefig(figure_path + '.pdf', bbox_inches='tight')

        return pivoted
