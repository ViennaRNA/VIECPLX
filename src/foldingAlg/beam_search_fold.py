#! /usr/bin/env python3

import RNA
import numpy as np
import math
import heapq
import argparse
import sys
import subprocess
from pathlib import Path

# from .basic import Folding
# from basic import Folding, Vienna
__all__ = ['BeamSearchFold']


#
# beam search algorithm below
#


class BeamSearch():

    # all potential base pairs (slow method)
    def canonical_pairs(self, seq):
        pairs = [[] for _ in seq]
        for i, nucleotide in enumerate(seq):
            for j in range(i+1, len(seq)):
                if (nucleotide == 'A' and seq[j] == 'U') or \
                    (nucleotide == 'U' and seq[j] == 'A') or \
                    (nucleotide == 'G' and seq[j] == 'C') or \
                    (nucleotide == 'C' and seq[j] == 'G') or \
                    (nucleotide == 'G' and seq[j] == 'U') or \
                        (nucleotide == 'U' and seq[j] == 'G'):
                    pairs[i].append(j+1)
        return pairs

    # potential base pairs using ViennaRNA base pair probabilities with defined threshold
    def pb_pairs(self, seq, threshold = 0.2):
        fc = RNA.fold_compound(seq)
        propensity, ensemble_energy = fc.pf()
        basepair_probs = fc.bpp()
        basepair_probs = np.array(basepair_probs)
        indices = np.where(basepair_probs > threshold)
        values = basepair_probs[indices]

        pairs = [[] for _ in seq]
        # print("Indices (i, j) and values above threshold:")
        for (i, j), value in zip(zip(*indices), values):
            # print(f"({i}, {j}) - {value}")
            pairs[int(i-1)].append(int(j))

        return pairs

    def predict_num_base_pairs(self, sequence_length, slope=0.306620, intercept=-3.622861, correction_factor=1.2):
        """
        Predicts the expected number of base pairs for a given RNA sequence length using a linear regression model.
        The prediction is based on linear regression of the form y = kx + d, 
        see beam_search_fold_playground.ipynb
        """
        predicted_base_pairs = slope * sequence_length + intercept
        predicted_base_pairs *= correction_factor
        # print (predicted_base_pairs)

        return math.ceil(predicted_base_pairs)



    def debug_print_structures(self, structures, rows=10, header=False):
        structure_length = structures[0][1][0]  # length in pt[0]

        energy_col_width = 7
        bp_col_width = 6
        struct_col_width = structure_length

        if type(rows) != int:
            rows = 10

        if header:
            header_line = f"{'bp':^{bp_col_width}} | {'Energy':^{energy_col_width}} | {'Structure':^{struct_col_width}}"
            separator = f"{'-' * bp_col_width} | {'-' * energy_col_width} | {'-' * struct_col_width}"
            print(header_line)
            print(separator)

        for energy, pt in structures[0:rows]:
            structure = RNA.db_from_ptable(pt)
            base_pairs = sum(1 for x in pt[1:] if x > 0) // 2
            print(
                f"{base_pairs:^{bp_col_width}} | {energy:^{energy_col_width}.2f} | {structure:<{struct_col_width}}")


    def fold(self, sequence, bp_threshold=0.01, bps=-1, l=40, debug=False):
        # Prepare the command to execute the C++ program with arguments
        command = [str(Path(Path(__file__).parent / 'cpp_experiments/beam_search_fold')), sequence, str(l), str(bp_threshold)]
        
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)       

        output = result.stdout.strip()
        return output



    def fold_legacy(self, seq, bp_threshold=0.01, bps=-1, l=10, debug=False):
        n = len(seq)
        s = '.' * n
        initial_pt = RNA.ptable(s)
        fc = RNA.fold_compound(seq)

        if bps==-1:
            bps = self.predict_num_base_pairs(n)

        if not bp_threshold:
            all_pairs = self.canonical_pairs(seq)
        else:
            all_pairs = self.pb_pairs(seq, threshold = bp_threshold)


        # Starting with open-chained structure
        structures = [(fc.eval_structure_pt(initial_pt), list(initial_pt))]
        best_structure = min(structures, key=lambda x: x[0])

        if debug:
            self.debug_print_structures(structures, rows=debug, header=True)

        for _ in range(bps):
            # new_structures = []
            # this is an ugly (Python) workaround to hash the pairing tables (converted into a tuple), 
            # such that only unique structures per iteration are kept. 
            new_structures_dict = {}
            for energy, pt in structures:
                # offspring structures per iteration
                l_indices = RNA.loopidx_from_ptable(pt)
                for i in range(1, n+1):
                    for j in all_pairs[i-1]:
                        if pt[i] == 0 and pt[j] == 0 and l_indices[i] == l_indices[j]:
                            en_candidate = fc.eval_move_pt(pt, i, j) + energy
                            # todo, pt_copy could be done after the sort to save some time
                            pt_copy = pt.copy()
                            pt_copy[i] = j
                            pt_copy[j] = i
                            # new_structures.append((en_candidate, pt_copy))
                            pt_tuple = tuple(pt_copy)
                            if pt_tuple not in new_structures_dict:
                                new_structures_dict[pt_tuple] = (en_candidate, pt_copy)


            new_structures = list(new_structures_dict.values())

            # sorting step
            # retain best l structures, sorted by energy (first tuple value)
            structures = heapq.nsmallest(l, new_structures, key=lambda x: x[0])
            if debug:
                self.debug_print_structures(
                    structures, rows=debug, header=False)

            if len(structures) == 0:
                break

            current_best = min(structures, key=lambda x: x[0])
            if current_best[0] < best_structure[0]:
                best_structure = current_best

        # return value is a dot-bracket string with energy in kcal/mol
        structure = RNA.db_from_ptable(best_structure[1])
        en = best_structure[0]/100.0

        # return structure, en
        return structure


def main():

    default_l = 20
    default_bps = -1

    parser = argparse.ArgumentParser(
        description='RNA Folding Simulation using Beam Search.')

    parser.add_argument('seq', type=str, help='Input RNA sequence to fold.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output.')
    parser.add_argument('-d', '--debug', default=False, help='Debug Settings')
    parser.add_argument('-l', type=int, default=default_l,
                        help=f'Beam search width (default: {default_l}).')
    parser.add_argument('-b', type=int, default=default_bps,
                        help=f'Max number of base pairs added (default: {default_bps}).')

    if len(sys.argv) == 1:
        parser.print_usage()
        print("\nSample call: python beam_search_fold.py AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGG -v -l 40")
        sys.exit(1)

    args = parser.parse_args()

    sequence = args.seq

    # BeamSearchFold call
    folding_algorithm = BeamSearch()
    # structure = folding_algorithm.fold(
    #     sequence, l=args.l, bps=args.b, debug=args.debug)
    structure = folding_algorithm.cpp_fold(
        sequence, l=args.l, bps=args.b)

    # vrna comparison
    
    s2 = RNA.fold(sequence)[0]

    fc = RNA.fold_compound(sequence)
    en = fc.eval_structure(structure)
    en2 = fc.eval_structure(s2)

    if args.verbose:
        print(f"RNA Sequence:     {sequence}")
        print(f"Folded Structure: {structure}{en:+6.2f} kcal/mol")
        print(f"Vienna RNA MFE:   {s2}{en2:+6.2f} kcal/mol")
    else:
        print(structure)


if __name__ == "__main__":
    main()
