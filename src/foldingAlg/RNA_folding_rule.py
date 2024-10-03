from RNA import eval_structure_simple
from RNA import random_string
from RNA import fold_compound
#import RNA



__all__ = ['RNAFoldingRule']

class Stem:
    def __init__(self, i, j, length,):
        self.i = i
        self.j = j
        self.length = length
        self.invalid = False
    def __str__(self):
        return f"Stem: i:{self.i}, j:{self.j}, l:{self.length}"



def fold(sequence, min_lenth = 1, min_dE = 0):
    """Folds the given RNA sequence similar to "An RNA Folding Rule" https://academic.oup.com/nar/article-abstract/12/1Part1/323/2889726

    Args:
        sequence (str): Sequence to fold
        min_lenth (int, optional): The minimal length of a helix for it to be considered when choosing which helices to create. Defaults to 3.
        min_dE (int, optional): Any helix hat would result in a higher delta energy when being formed is not considered for folding. Defaults to 0.

    Returns:
        str: DotBracket of the resulting secondary structure
    """

    fc = fold_compound(sequence)
    current_ener = 0
    helices_per_nuc = []
    db = '.'*len(sequence)
    for _ in sequence:
        helices_per_nuc.append([])

    helices = find_all_helices(fc,sequence,min_lenth,helices_per_nuc,db)
    helices = dict(sorted(helices.items(), key=lambda helix: helix[1]))




    #Insert helices until no more insertions are possible
    while len(helices) != 0:
        #for s in helices:
            #print(str(s) + ": " + str(helices[s]))
        #print("\n")

        current_helix = next(iter(helices))
        if helices[current_helix] > min_dE: break
        db = apply_to_db(db,current_helix)
        remove_crossing(helices,current_helix,helices_per_nuc)
        update_ener(fc,helices,db,sequence,current_helix)
        helices = dict(sorted(helices.items(), key=lambda helix: helix[1]))


    return ''.join(db)





def find_all_helices(fc, sequence, min_lenth, helices_per_nuc, db):
    helices = {}
    for j in range(len(sequence)):
        for i in range(j):

            helix_length = 0

            while bases_compatible(sequence,i+helix_length,j-helix_length):
                helix_length += 1

            if helix_length < min_lenth:
               continue

            new_helix = Stem(i,j,helix_length)
            helices[new_helix] = get_ener(fc,new_helix,db,sequence)

            for i in range(new_helix.i,new_helix.i+new_helix.length):
                helices_per_nuc[i].append(new_helix)

            for i in range(new_helix.j,new_helix.j-new_helix.length,-1):
                helices_per_nuc[i].append(new_helix)

    return helices

def remove_crossing(helixs,current_helix,helices_per_nuc):

    #remove all helices that have overlapping bases with the new helix
    for i in range(current_helix.i,current_helix.i+current_helix.length):
        for possible_crossing in helices_per_nuc[i]:
            if possible_crossing.invalid: continue 
            del helixs[possible_crossing]
            possible_crossing.invalid = True
        helices_per_nuc[i].clear()

    #remove all helices that have overlapping bases with the new helix
    for i in range(current_helix.j,current_helix.j-current_helix.length,-1):
        for possible_crossing in helices_per_nuc[i]:
            if possible_crossing.invalid: continue 
            del helixs[possible_crossing]
            possible_crossing.invalid = True
        helices_per_nuc[i].clear()
    
    #remove bases that have bases within and outside of the new helix
    for i in range(current_helix.i+current_helix.length,(current_helix.j-current_helix.length)+1):
        index_for_nuc = 0
        while index_for_nuc <= len(helices_per_nuc[i])-1:
        #for possible_crossing in helices_per_nuc[i]:
            possible_crossing = helices_per_nuc[i][index_for_nuc]
            if possible_crossing.invalid:
                helices_per_nuc[i] = helices_per_nuc[i][:index_for_nuc] + helices_per_nuc[i][index_for_nuc+1 :]
                continue 
            
            length_until = min(i-possible_crossing.i,possible_crossing.j-i)
            if pos_in_loop(current_helix,possible_crossing.i+length_until) and pos_in_loop(current_helix,possible_crossing.j-length_until):
                index_for_nuc += 1
                continue 

            del helixs[possible_crossing]
            possible_crossing.invalid = True
            helices_per_nuc[i] = helices_per_nuc[i][:index_for_nuc] + helices_per_nuc[i][index_for_nuc+1 :]

def update_ener(fc,helices,db,seq,current_helix):
    for helix in helices:
        helices[helix] = get_ener(fc,helix,db,seq)

def get_ener(fc,helix,db,seq):
    dE = 0
    db_cp = db
    for curr_length in range(0,helix.length):
        bp_i = helix.i+curr_length
        bp_j = helix.j-curr_length


        dE += fc.eval_move(db_cp,bp_i+1,bp_j+1)
        db_cp = db_cp[:bp_i] + '(' + db_cp[bp_i + 1:bp_j] + ')' + db_cp[bp_j+1:]
    return dE

def pos_in_loop(helix,i):
    if i >= helix.i+helix.length and i <= helix.j-helix.length:
        return True
    return False


def bases_compatible(seq, base_1_index, base_2_index):
    if base_1_index < 0 or base_1_index >= len(seq):
        return False
    if base_2_index < 0 or base_2_index >= len(seq):
        return False
    base_1 = seq[base_1_index]
    base_2 = seq[base_2_index]
    
    if abs(base_2_index - base_1_index) <= 3:
        return False

    if base_1 in {'A','a'} and base_2 in {'U', 'u'}:
        return True
    if base_1 in {'U','u'} and base_2 in {'A', 'a'}:
        return True
    if base_1 in {'C','c'} and base_2 in {'G', 'g'}:
        return True
    if base_1 in {'G','g'} and base_2 in {'C', 'c'}:
        return True
    if base_1 in {'U','u'} and base_2 in {'G', 'g'}:
        return True
    if base_1 in {'G','g'} and base_2 in {'U', 'u'}:
        return True
    return False

def apply_to_db(db,helix):
    for curr_length in range(0,helix.length):
        bp_i = helix.i+curr_length
        bp_j = helix.j-curr_length

        assert(db[bp_i] == '.' and db[bp_j] == '.')

        db = db[:bp_i] + '(' + db[bp_i + 1:bp_j] + ')' + db[bp_j+1:]
    return db



def main():
      
    seq = random_string(100,"ACGU")
    db = fold(seq)
    print(seq + " " + db + " " + str(eval_structure_simple(seq,db,0,None)))

#   nr = 10000
#   for i in range(1,10):
#      sum_distance = 0
#      outF = [RNA.random_string(100,"ACGU") for i in range(nr)]
#      for seq in outF:
#          db = fold(seq,i,0)

#          fc  = RNA.fold_compound(seq)
#          (ss, mfe) = fc.mfe()
#          algo_tree = RNA.db_to_tree_string(db,RNA.STRUCTURE_TREE_SHAPIRO)
#          RNAfold_tree = RNA.db_to_tree_string(ss,RNA.STRUCTURE_TREE_SHAPIRO)
#          algo_tree = RNA.make_tree(algo_tree)
#          RNAfold_tree = RNA.make_tree(RNAfold_tree)

#          sum_distance += RNA.tree_edit_distance(algo_tree,RNAfold_tree)
#          #print(seq + " " + db + " " + str(eval_structure_simple(seq,db,0,None)))
#      print(f"{i}: {sum_distance/nr}")


if __name__ == "__main__":
    main()
 
