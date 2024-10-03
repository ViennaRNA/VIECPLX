import random
from itertools import islice

# Define complement nucleotides and translation for folding representation.
# Added optional nucleotides for G and U to allow for more complex folding.
complements = {'A': 'U', 'U': 'AG', 'C': 'G', 'G': 'CU'}
translate = {0: '.', 1: '(', 2: ')'}
min_loop_size = 3

DEBUG = False

def fold(rna_sequence: str, window_size: int = 3) -> str:
    """
    Simulate the folding of an RNA sequence.

    Args:
    rna_sequence (str): The RNA sequence to be folded.
    window_size (int): The size of the window to look for complement nucleotides. Defaults to 3.

    Returns:
    str: A string representing the folded RNA, where '.' indicates an unpaired nucleotide, 
    '(' a nucleotide paired at the start, and ')' a nucleotide paired at the end.
    """
    if not rna_sequence:
        return ''
    stack = [(0, len(rna_sequence))]
    folded_string = [0 for _ in rna_sequence]
    while stack:
        start_index, end_index = stack.pop()
        find_subset(rna_sequence, start_index, end_index, window_size, stack, folded_string)
    return ''.join(translate[n] for n in folded_string)

def find_subset(rna_seq: str, start_index: int, end_index: int, window_size: int, stack: list, folded_string: list):
    """
    Find a subset of the RNA sequence that can be folded and update the folding status.

    Args:
    rna_seq (str): The RNA sequence.
    start_index (int): The starting index of the current segment.
    end_index (int): The ending index of the current segment.
    window_size (int): The window size for folding.
    stack (list): The stack to manage different segments.
    folded_string (list): The list representing the folding status of the RNA sequence.
    """
    for index in range(start_index, end_index - window_size + 1):
        window = rna_seq[index:index + window_size]


        reverse_search_start_index = max(0, index - min_loop_size)
        search_range = rna_seq[start_index:reverse_search_start_index]

        reverse_complement_result = reverse_complement(window)

        replacement_index = None

        for (window_index, window_slice) in enumerate(windows(search_range, window_size)):
            if all([base in bases for (base, bases) in zip(window_slice,reverse_complement_result)]):
                replacement_index = start_index + window_index
                break
                
        if replacement_index != None:
            replace_in_range(folded_string, 1, replacement_index, replacement_index + window_size)

            replace_in_range(folded_string, 2, index, index + window_size)

            zip_option = zip_nucleotides_from_seed(rna_seq,  replacement_index, index, window_size, folded_string)
            if zip_option:
                stack.append((zip_option[1], end_index))
            else:
                stack.append((index + window_size, end_index))
            return

def zip_nucleotides_from_seed(rna_seq: str, seed_start: int, seed_end: int, seed_length: int, folded_string: list) -> tuple:
    """
    Zip nucleotides outwards and inwards from a seed region in the RNA sequence.
    Args:
    rna_seq (str): The RNA sequence.
    seed_start (int): The start index of the seed region.
    seed_end (int): The end index of the region to zip with the seed.
    seed_length (int): The length of the seed region.
    folded_string (list): The list representing the folding status of the RNA sequence.
    Returns:
    tuple: A tuple indicating the next segment to be analyzed or None if no further zipping is possible.
    """
    pointer_low_in = seed_start + seed_length
    pointer_high_in = seed_end - 1

    pointer_low_out = seed_start - 1
    pointer_high_out = seed_end + seed_length
    
    if DEBUG:
        print(folded_string)
        print(f"Inward pointers: {pointer_low_in}, {pointer_high_in}")
        print(f"Outward pointers: {pointer_low_out}, {pointer_high_out}")
    
    last_zip_pointer_pair = None
    
    stack_size = 0
    # Zip inwards
    while True:

        zipped = False
        # This is probably not necessary
        if pointer_low_in < pointer_high_in - min_loop_size:
            if rna_seq[pointer_low_in] in complements[rna_seq[pointer_high_in]]:
                stack_size += 1
                if stack_size == 2:
                    replace_in_range(folded_string, 1, pointer_low_in, pointer_low_in + 1)
                    replace_in_range(folded_string, 2, pointer_high_in - 1, pointer_high_in)
                elif stack_size > 2:
                    replace_index(folded_string, pointer_low_in, 1)
                    replace_index(folded_string, pointer_high_in, 2)
                else:
                    stack_size = 0

            pointer_low_in += 1
            pointer_high_in -= 1
        else:
            break

    while True:
        zipped = False
        
        # Zip outwards
        if pointer_low_out >= 0 and pointer_high_out < len(rna_seq):
            if rna_seq[pointer_low_out] in complements[rna_seq[pointer_high_out]]:
                if folded_string[pointer_low_out] == 0 and folded_string[pointer_high_out] == 0:
                    replace_index(folded_string, pointer_low_out, 1)
                    replace_index(folded_string, pointer_high_out, 2)
                    zipped = True
            pointer_low_out -= 1
            pointer_high_out += 1
        
        if not zipped:
            last_zip_pointer_pair = (pointer_low_out, pointer_high_out)
            break
    
    return last_zip_pointer_pair if last_zip_pointer_pair else None

def replace_in_range(string: list, character: int, start: int, end: int):
    """
    Replace characters in a range of a list with a given character.

    Args:
    string (list): The list to be modified.
    character (int): The character to replace with.
    start (int): The starting index of the range.
    end (int): The ending index of the range.
    """
    for index in range(start, end):
        if string[index] == 0:
            string[index] = character


def replace_index(string: list, index: int, character: int):
    if string[index] == 0:
        string[index] = character

def windows(seq, n=3):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def reverse_complement(seq: str) -> str:
    """
    Get the reverse complement of an RNA sequence.

    Args:
    seq (str): The RNA sequence.

    Returns:
    str: The reverse complement of the given RNA sequence.
    """
    return [complements[base] for base in reversed(seq)]
