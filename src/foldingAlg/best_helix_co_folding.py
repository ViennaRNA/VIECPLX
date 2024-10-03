from RNA import random_string  # Only used to generate random sequences
from RNA import eval_structure_pt_simple
from RNA import eval_structure_simple
from RNA import fold_compound
from RNA import CONSTRAINT_DB_DEFAULT
import RNA


__all__ = ['BestHelixCoFolding']


def fold(seq: str, helix_length: int = 6, helix_length2: int = 4, min_loop_size: int = 3, helix_better_than: float = 0, growth: int = 3, useRNAfoldForRound2=False):
    """Fold the given RNA Sequence and returns a dotbracket. In each step, grow the chain and then form the best helix that can be
    constructed using the newly inserted base. This procedure is done two times with a different minimum helix length. Alternatively,
    the second stage can be replaced by running the RNAFold algorithm with the structure from the first stage as constraints.

    Args:
        seq (str): The RNA sequence in Upper case letters
        min_helix_length1 (int, optional): The minimal length a helix must have in order to be formed in stage 1. Defaults to 3.
        min_helix_length2 (int, optional): The minimal length a helix must have in order to be formed in stage 2, if RNAFold is not used. Defaults to 3.
        min_loop_size (int, optional): The minimal number of unpaired bases between two paired bases. Defaults to 3.
        helix_better_than (float, optional): Only use helices which result in a structure with lower energy than this. Defaults to 0.
        growth (int, optional): Each step, grow the chain by this value and then seach for the best helix which contains at least one of the new bases. Defaults to 3.
        useRNAfoldForRound2 (bool, optional): Substitute the second stage with running the RNAFold algorithm using the structure of the first stage as constraints. Defaults to False

    Returns:
        (str): Dotbracket sting of the calculated fold
    """
    if not seq:
        return ""
    seq_length = len(seq)
    seq = ' ' + seq
    pt = [0] * (seq_length+1)
    pt[0] = seq_length

    pt = iterate_i(seq, pt, helix_length, min_loop_size, helix_better_than, comparator_energy, growth, seq_length)
    # return pt_to_db(pt)

    if useRNAfoldForRound2:
        fc = fold_compound(seq[1:])
        # print(pt_to_db(pt))
        fc.constraints_add(pt_to_db(pt), CONSTRAINT_DB_DEFAULT)
        (ss, mfe) = fc.mfe()
        return ss

    else:
        pt = iterate_i(seq, pt, helix_length2, min_loop_size, helix_better_than, comparator_energy, growth, seq_length)
        return pt_to_db(pt)


def iterate_i(seq: str, pt, min_helix_length: int, min_loop_size: int, helix_better_than: int, comperator, growth: int, seq_length: int):
    position_i = min(growth, seq_length)
    while position_i <= seq_length:
        best_helix_start, best_helix_end, used_i = return_best_helix(seq, pt, position_i, min_helix_length, min_loop_size, helix_better_than, growth, comperator)
        if best_helix_start is not None:
            pt = create_helix(pt, best_helix_start, best_helix_end, used_i)

        if (position_i != seq_length) and (position_i + growth > seq_length):
            growth = seq_length - position_i
            position_i = seq_length
        else:
            position_i += growth

    return pt


def return_best_helix(seq: str, pt, position_i: int, min_helix_length: int, min_loop_size: int, helix_better_than: float, growth: int, comparator):
    best_helix_start = None
    best_helix_end = None
    best_helix_value = None
    best_i = None
    helix_start = 0
    helix_end = 0
    current_i = position_i

    for i in range(0, growth):
        position_j = position_i-1-i
        if pt[current_i] != 0:
            current_i -= 1
            position_j -= 1
            continue

        while position_j > 0:
            if pt[position_j] < position_j and pt[position_j] != 0:
                position_j = pt[position_j]-1
                continue
            if pt[position_j] > current_i:
                break
            if not bases_compatible(seq, position_j, current_i, min_loop_size):
                position_j -= 1
                continue

            assert (pt[position_j] == 0)

            helix_start = position_j
            helix_end = position_j+1

            while helix_end < current_i-(helix_end-helix_start):
                current_length = (helix_end-helix_start)
                if pt[helix_end] != 0 or pt[current_i-current_length] != 0:
                    break
                if not bases_compatible(seq, helix_end, current_i-current_length, min_loop_size):
                    break

                helix_end += 1

            helix_length = helix_end-helix_start
            if helix_length < min_helix_length:
                position_j -= 1
                continue
            is_better, value = comparator(seq, pt, helix_start, helix_end, current_i, best_helix_value, helix_better_than)
            if is_better:
                best_helix_start = helix_start
                best_helix_end = helix_end
                best_i = current_i
                best_helix_value = value
            position_j -= 1

        current_i -= 1
        position_j -= 1

    return best_helix_start, best_helix_end, best_i


def comparator_length(seq: str, pt, helix_start: int, helix_end: int, position_i: int, best_helix_value, helix_better_than: float):
    value = helix_end - helix_start
    if helix_better_than is not None and value > helix_better_than:
        return False, None
    if best_helix_value is None:
        return True, value
    if value > best_helix_value:
        return True, value
    return False, value


def comparator_energy(seq: str, pt, helix_start: int, helix_end: int, position_i: int, best_helix_value, helix_better_than: float):
    value = eval_structure_pt_simple(seq[1:], create_helix(pt, helix_start, helix_end, position_i), 0, None) / 100
    if helix_better_than is not None and value > helix_better_than:
        return False, None
    if best_helix_value is None:
        return True, value
    if value < best_helix_value:
        return True, value
    return False, value


def create_helix(pt, best_helix_start: int, best_helix_end: int, position_i: int):
    pt_new = pt.copy()
    for j in range(best_helix_start, best_helix_end):
        pt_new[j] = position_i - (j - best_helix_start)
        pt_new[position_i - (j - best_helix_start)] = j
    return pt_new


def bases_compatible(seq: str, base_1_index: int, base_2_index: int, min_loop_size: int):
    base_1 = seq[base_1_index]
    base_2 = seq[base_2_index]

    if abs(base_2_index - base_1_index) <= min_loop_size:
        return False

    if base_1 == 'A' and base_2 == 'U':
        return True
    if base_1 == 'U' and base_2 == 'A':
        return True
    if base_1 == 'C' and base_2 == 'G':
        return True
    if base_1 == 'G' and base_2 == 'C':
        return True
    if base_1 == 'U' and base_2 == 'G':
        return True
    if base_1 == 'G' and base_2 == 'U':
        return True
    return False


def pt_to_db(pt):
    db = ""
    for i in range(1, pt[0]+1):
        if pt[i] == 0:
            db += '.'
        elif pt[i] > i:
            db += '('
        elif pt[i] < i:
            db += ')'
    return db


def check_if_pt_valid(pt):
    for base_i, paired_i in enumerate(pt[1:]):
        base_i += 1
        if paired_i < base_i:
            continue
        for base_j in range(base_i, paired_i):
            if pt[base_j] == 0:
                continue
            if pt[base_j] > paired_i or pt[base_j] < base_i:
                return False
    return True


def main():

    for _ in range(20):
        seq = random_string(100, "ACGU")
        # seq = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU"
        db = fold(seq)
        print(seq + " " + db + " " + str(eval_structure_simple(seq, db, 0, None)))

#   nr = 500

#   for i in range(1, 10):
#       for r in range(1, i+1):
#           sum_distance = 0
#           outF = [RNA.random_string(100, "ACGU") for i in range(nr)]
#           j = 0
#           for seq in outF:
#               if j == nr:
#                   break
#               db = fold(seq, i, r, useRNAfoldForRound2=False)
#               j += 1

#               fc = RNA.fold_compound(seq)
#               (ss, mfe) = fc.mfe()
#               algo_tree = RNA.db_to_tree_string(db, RNA.STRUCTURE_TREE_SHAPIRO)
#               RNAfold_tree = RNA.db_to_tree_string(ss, RNA.STRUCTURE_TREE_SHAPIRO)
#               algo_tree = RNA.make_tree(algo_tree)
#               RNAfold_tree = RNA.make_tree(RNAfold_tree)

#               sum_distance += RNA.tree_edit_distance(algo_tree, RNAfold_tree)
#               # print(seq + " " + db + " " + str(eval_structure_simple(seq,db,0,None)))
#           print(f"{i},{r}: {sum_distance/nr}")


if __name__ == "__main__":
    main()
