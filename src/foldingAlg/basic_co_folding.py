from RNA import random_string  # Only used to generate random sequences
from RNA import eval_structure_simple
from RNA import fold_compound
from RNA import CONSTRAINT_DB_DEFAULT
# import RNA


def fold(seq: str, helix_length: int = 6, helix_length2: int = 4, direction_first_base_i: str = "left_to_right", direction_second_base_j: str = "left_to_i", use_long_helix: bool = True, useRNAfoldForRound2=False):
    """Folds the given RNA-Sequence

    Args:
        seq (str): Sequence of the RNA to fold in capital letters
        helix_length (int, optional): fixed length of the created helices. Defaults to 1.
        direction_first_base_i (str, optional): How to loop over the bases to find the first base of the pair (i,j). Can be "left_to_right" or "right_to_left". Defaults to "left_to_right".
        direction_second_base_j (str, optional): How to loop over the bases to find the second base of the pair (i,j). Can be "<x>_to_<y>", where x and y can be "left", "right", "i"; "i" is the base already chosen in the first loop. Defaults to "left_to_i".
        use_long_helix (bool, optional): If this is set to true, the first helix with as many bases as possible but at least as many bases as speciefied in "helix_length" will be paired. If set to false, only helices with the exact length as defined by the arguments of this function will be formed. Defaults to True.



    Returns:
        str: 2D structure in DotBracket notation
    """
    if not seq:
        return ""

    seq_length = len(seq)
    seq = ' ' + seq
    pt = [0] * (seq_length+1)
    pt[0] = seq_length

    iterate_i_func = None
    iterate_j_func = None
    stop_at_i = False
    start_at_i = False

    match direction_first_base_i:
        case "left_to_right":
            iterate_i_func = iterate_over_i_left_to_right
        case "right_to_left":
            iterate_i_func = iterate_over_i_right_to_left
        case _:
            raise Exception("Invalid opion for 'direction_first_base_i'")

    match direction_second_base_j:
        case "left_to_i":
            iterate_j_func = iterate_over_j_left_to_right
            stop_at_i = True
            start_at_i = False
        case "left_to_right":
            iterate_j_func = iterate_over_j_left_to_right
            stop_at_i = False
            start_at_i = False
        case "right_to_i":
            iterate_j_func = iterate_over_j_right_to_left
            stop_at_i = True
            start_at_i = False
        case "right_to_left":
            iterate_j_func = iterate_over_j_right_to_left
            stop_at_i = False
            start_at_i = False
        case "i_to_right":
            iterate_j_func = iterate_over_j_left_to_right
            stop_at_i = False
            start_at_i = True
        case "i_to_left":
            iterate_j_func = iterate_over_j_right_to_left
            stop_at_i = False
            start_at_i = True
        case _:
            raise Exception("Invalid opion for 'direction_second_base_j'")

    iterate_i_func(seq_length, helix_length, seq, pt, iterate_j_func, stop_at_i, start_at_i, use_long_helix)

    if useRNAfoldForRound2:
        fc = fold_compound(seq[1:])
        # print(pt_to_db(pt))
        fc.constraints_add(pt_to_db(pt), CONSTRAINT_DB_DEFAULT)
        (ss, mfe) = fc.mfe()
        return ss

    else:
        iterate_i_func(seq_length, helix_length2, seq, pt, iterate_j_func, stop_at_i, start_at_i, use_long_helix)
    db = pt_to_db(pt)
    return db


def iterate_over_j_left_to_right(seq_length: int, helix_length: int, helix_direction_i: int, seq: str, pt, i: int, stop_at_i: bool, start_at_i: bool, use_long_helix: bool):
    """Iterate over the sequence from left to right to search for a base j, which can pair with base i

    Args:
        seq_length (int): Length of the complete Sequence
        helix_length (int): Exact length of the helix we want to create
        helix_direction_i (int): Should be 1 if the bases from i to i+helix_lenght should be paired and -1 if the bases from i-helix_lengh to i should be paired
        seq (str): Complete sequence
        pt (_type_): Pairtable of the RNA
        stop_at_i (bool): True if the searches for j should stop if it passes over base i
        start_at_i (bool): True if the search for j should start at the position of base i
        use_long_helix (bool, optional): If this is set to true, the first helix with as many bases as possible but at least as many bases as speciefied in "helix_length" will be paired. If set to false, only helices with the exact length as defined by the arguments of this function will be formed.

    Returns:
        bool: True if a helix was created and False otherwise
    """
    if pt[i] != 0:
        return False
    j: int = i

    # If we should start from the left side, move from position i to the left until
    # we are either on the far left or we encounter a basepair that crosses i. In
    # the second case, only consider bases within the crossing bp to avoid pseudoknots
    if not start_at_i:
        while j > 1 and pt[j-1] < i:
            j -= 1
            if j == 1:
                break

    while j <= seq_length-(helix_length-1):
        if stop_at_i and j >= i-(helix_length-1):
            break
        if (not use_long_helix) and pair_if_possible(seq_length, helix_length, helix_direction_i, 1, seq, pt, i, j):
            return True
        if use_long_helix and pair_longest_helix_possible(seq_length, helix_length, helix_direction_i, -helix_direction_i, seq, pt, i, j):
            return True
        if pt[j] > j:
            j = pt[j]+1
        elif pt[j] == 0:
            j += 1
        else:
            break
    return False


def iterate_over_j_right_to_left(seq_length: int, helix_length: int, helix_direction_i: int, seq: str, pt, i: int, stop_at_i: bool, start_at_i: bool, use_long_helix: bool):
    """Iterate over the sequence from right to left to search for a base j, which can pair with base i

    Args:
        seq_length (int): Length of the complete Sequence
        helix_length (int): Exact length of the helix we want to create
        helix_direction_i (int): Should be 1 if the bases from i to i+helix_lenght should be paired and -1 if the bases from i-helix_lengh to i should be paired
        seq (str): Complete sequence
        pt (_type_): Pairtable of the RNA
        stop_at_i (bool): True if the searches for j should stop if it passes over base i
        start_at_i (bool): True if the search for j should start at the position of base i
        use_long_helix (bool, optional): If this is set to true, the first helix with as many bases as possible but at least as many bases as speciefied in "helix_length" will be paired. If set to false, only helices with the exact length as defined by the arguments of this function will be formed.

    Returns:
        bool: True if a helix was created and False otherwise
    """
    if pt[i] != 0:
        return False
    j: int = i

    # If we should start from the right side, move from position i to the right until
    # we are either on the far right or we encounter a basepair that crosses i. In
    # the second case, only consider bases within the crossing bp to avoid pseudoknots
    if not start_at_i:
        while j < seq_length and (pt[j+1] > i or pt[j+1] == 0):
            j += 1
            if j == seq_length:
                break

    while j >= 1+(helix_length-1):
        if stop_at_i and j <= i+(helix_length-1):
            break
        if (not use_long_helix) and pair_if_possible(seq_length, helix_length, helix_direction_i, -1, seq, pt, i, j):
            return True
        if use_long_helix and pair_longest_helix_possible(seq_length, helix_length, helix_direction_i, -helix_direction_i, seq, pt, i, j):
            return True

        if pt[j] < j and pt[j] != 0:
            j = pt[j]-1
        elif pt[j] > j:
            break
        else:
            j -= 1
    return False


def iterate_over_i_left_to_right(seq_length: int, helix_length: int, seq: str, pt, iterate_over_j, stop_at_i: bool, start_at_i: bool, use_long_helix: bool):
    """Iterate over all bases from left to right. Call this base i. Call a function that searches for a base j, which can be paired with i

    Args:
        seq_length (int): Length of the complete sequence
        helix_length (int): Exact length of the helix we want to create
        seq (str): Complete sequence
        pt (_type_): Pairtable of the RNA
        iterate_over_j (_type_): Function that searches for the second base j to pair with base i
        stop_at_i (bool): True if the function that searches for j should stop if it passes over base i
        start_at_i (bool): True if the function that searches for j should start at the position of base i
        use_long_helix (bool, optional): If this is set to true, the first helix with as many bases as possible but at least as many bases as speciefied in "helix_length" will be paired. If set to false, only helices with the exact length as defined by the arguments of this function will be formed.

    """
    i: int = 1
    while i <= seq_length-(helix_length-1):
        if pt[i] == 0:
            iterate_over_j(seq_length, helix_length, 1, seq, pt, i, stop_at_i, start_at_i, use_long_helix)
        i += 1


def iterate_over_i_right_to_left(seq_length: int, helix_length: int, seq: str, pt, iterate_over_j, stop_at_i: bool, start_at_i: bool, use_long_helix: bool):
    """Iterate over all bases from right to left. Call this base i. Call a function that searches for a base j, which can be paired with i

    Args:
        seq_length (int): Length of the complete sequence
        helix_length (int): Exact length of the helix we want to create
        seq (str): Complete sequence
        pt (_type_): Pairtable of the RNA
        iterate_over_j (_type_): Function that searches for the second base j to pair with base i
        stop_at_i (bool): True if the function that searches for j should stop if it passes over base i
        start_at_i (bool): True if the function that searches for j should start at the position of base i
        use_long_helix (bool, optional): If this is set to true, the first helix with as many bases as possible but at least as many bases as speciefied in "helix_length" will be paired. If set to false, only helices with the exact length as defined by the arguments of this function will be formed.
    """
    i: int = seq_length
    while i >= 1+(helix_length-1):
        if pt[i] == 0:
            iterate_over_j(seq_length, helix_length, -1, seq, pt, i, stop_at_i, start_at_i, use_long_helix)
        i -= 1


def pair_if_possible(seq_length: int, helix_length: int, helix_direction_i: int, helix_direction_j: int, seq: str, pt, i: int, j: int):
    """Pairs bases to a helix of a given size if possible

    Args:
        seq_length (int): Length of the complete sequence
        helix_length (int): Exact length of the helix that can be paired
        helix_direction_i (int): Should be 1 if the bases from i to i+helix_lenght should be paired and -1 if the bases from i-helix_lengh to i should be paired
        helix_direction_j (int): Should be 1 if the bases from j to j+helix_lenght should be paired and -1 if the bases from j-helix_lengh to j should be paired
        seq (str): Complete Sequence of the RNA
        pt (_type_): Pairtable of the complete RNA
        i (int): beginning of the first segment that should be paired
        j (int): beginning of the second segment that should be paired

    Returns:
        bool: True if it was possible to create a helix and False otherwise
    """
    if interval_compatible(seq, pt, i, j, helix_length, helix_direction_i, helix_direction_j):
        pair_helix(pt, i, j, helix_length, helix_direction_i, helix_direction_j)
        return True
    return False


def pair_longest_helix_possible(seq_length: int, min_helix_length: int, helix_direction_i: int, helix_direction_j: int, seq: str, pt, i: int, j: int):
    """Pairs bases to a helix of a given size if possible

    Args:
        seq_length (int): Length of the complete sequence
        min_helix_length (int): Minimal length of the helix that can be paired
        helix_direction_i (int): Should be 1 if the bases from i to i+helix_lenght should be paired and -1 if the bases from i-helix_lengh to i should be paired
        helix_direction_j (int): Should be 1 if the bases from j to j+helix_lenght should be paired and -1 if the bases from j-helix_lengh to j should be paired
        seq (str): Complete Sequence of the RNA
        pt (_type_): Pairtable of the complete RNA
        i (int): beginning of the first segment that should be paired
        j (int): beginning of the second segment that should be paired

    Returns:
        bool: True if it was possible to create a helix and False otherwise
    """
    # Try in which direction we can generate a longer helix and choose the longest possbile helix
    possible_helix_length_first_direction = longest_interval_compatible(seq, pt, i, j, helix_direction_i, helix_direction_j)
    possible_helix_length_second_direction = longest_interval_compatible(seq, pt, i, j, -helix_direction_i, -helix_direction_j)
    if possible_helix_length_first_direction + possible_helix_length_second_direction > min_helix_length:
        pair_helix(pt, i, j, possible_helix_length_first_direction, helix_direction_i, helix_direction_j)
        pair_helix(pt, i, j, possible_helix_length_second_direction, -helix_direction_i, -helix_direction_j)
        return True
    # if possible_helix_length_first_direction >= min_helix_length and possible_helix_length_first_direction > possible_helix_length_second_direction:
        # pair_helix(pt,i,j,possible_helix_length_first_direction,helix_direction_i,helix_direction_j)
        # return True
    # if possible_helix_length_second_direction >= min_helix_length and possible_helix_length_second_direction > possible_helix_length_first_direction:
        # pair_helix(pt,i,j,possible_helix_length_second_direction,-helix_direction_i,-helix_direction_j)
        # return True
    return False


def pair_helix(pt, i: int, j: int, helix_length: int, helix_direction_i: int, helix_direction_j: int):
    for offset in range(0, helix_length):
        pt[i+offset*helix_direction_i] = j+offset*helix_direction_j
        pt[j+offset*helix_direction_j] = i+offset*helix_direction_i


def interval_compatible(seq: str, pt, i: int, j: int, helix_length: int, helix_direction_i: int, helix_direction_j: int):
    for offset in range(0, helix_length):
        if not bases_compatible(seq, i+offset*helix_direction_i, j+offset*helix_direction_j):
            return False
        if pt[i+offset*helix_direction_i] != 0 or pt[j+offset*helix_direction_j] != 0:
            return False
    return True


def longest_interval_compatible(seq: str, pt, i: int, j: int, helix_direction_i: int, helix_direction_j: int):
    offset = 0
    while bases_compatible(seq, i+offset*helix_direction_i, j+offset*helix_direction_j) and pt[i+offset*helix_direction_i] == 0 and pt[j+offset*helix_direction_j] == 0:
        offset += 1

    return offset


def bases_compatible(seq: str, base_1_index: int, base_2_index: int):
    if base_1_index < 1 or base_1_index >= len(seq):
        return False
    if base_2_index < 1 or base_2_index >= len(seq):
        return False
    base_1 = seq[base_1_index]
    base_2 = seq[base_2_index]

    if abs(base_2_index - base_1_index) <= 3:
        return False

    if base_1 in {'A', 'a'} and base_2 in {'U', 'u'}:
        return True
    if base_1 in {'U', 'u'} and base_2 in {'A', 'a'}:
        return True
    if base_1 in {'C', 'c'} and base_2 in {'G', 'g'}:
        return True
    if base_1 in {'G', 'g'} and base_2 in {'C', 'c'}:
        return True
    if base_1 in {'U', 'u'} and base_2 in {'G', 'g'}:
        return True
    if base_1 in {'G', 'g'} and base_2 in {'U', 'u'}:
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

    seq = random_string(100, "ACGU")
    seq = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU"
    db = fold(seq, useRNAfoldForRound2=True)
    print(seq + " " + db + " " + str(eval_structure_simple(seq, db, 0, None)))

#   nr = 1000
#   for i in range(3,10):
#       for r in range(1,2):
#          sum_distance = 0
#          outF = [RNA.random_string(100,"ACGU") for i in range(nr)]
#          j = 0
#          for seq in outF:
#              if j == nr: break
#              db = fold(seq,i,useRNAfoldForRound2=True)
#              j += 1

#              fc  = RNA.fold_compound(seq)
#              (ss, mfe) = fc.mfe()
#              algo_tree = RNA.db_to_tree_string(db,RNA.STRUCTURE_TREE_SHAPIRO)
#              RNAfold_tree = RNA.db_to_tree_string(ss,RNA.STRUCTURE_TREE_SHAPIRO)
#              algo_tree = RNA.make_tree(algo_tree)
#              RNAfold_tree = RNA.make_tree(RNAfold_tree)

#              sum_distance += RNA.tree_edit_distance(algo_tree,RNAfold_tree)
#              #print(seq + " " + db + " " + str(eval_structure_simple(seq,db,0,None)))
#          print(f"{i},{r}: {sum_distance/nr}")


if __name__ == "__main__":
    main()
