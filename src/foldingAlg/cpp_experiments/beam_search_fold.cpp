/*
Compilation:

Option 1: If ViennaRNA dev header files are present in the system's default path (e.g., /usr/include/)
g++ beam_search_fold.cpp -o beam_search_fold -lRNA -lm -lstdc++ -fopenmp -std=c++20

Option 2: If ViennaRNA is installed inside a Conda environment (viennarna package via bioconda)
g++ beam_search_fold.cpp -o beam_search_fold \
-I/scr/XXX/YYY/miniforge/envs/YOURCONDAENV/include \
-L/scr/XXX/YYY/miniforge/envs/YOURCONDAENV/lib \
-lRNA -lm -lstdc++ -fopenmp -std=c++20

After compiling, change permissions:
chmod +x beam_search_fold
*/



#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <set>
#include <map>
#include <cstring>
#include <random>
#include <utility>
#include <unordered_set>

// 3rd party header only libraries
#include "libraries/robin_set.h"
#include "libraries/robin_map.h"

extern "C"
{
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/landscape/findpath.h"
#include "ViennaRNA/landscape/neighbor.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/structures.h"
#include <ViennaRNA/part_func.h>
}

void debug_print_moves(const std::vector<std::pair<int, int>> &vec)
{
    std::cout << "[";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        std::cout << "(" << vec[i].first << ", " << vec[i].second << ")";
        if (i < vec.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

uint64_t int_hash_64(uint64_t x)
{
    uint64_t z = (x += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

struct intermediate
{
    std::vector<std::pair<int, int>> moves;
    int curr_en;
    uint64_t s_hash;

    // Default constructor (only needed for vector of intermediate resize)
    intermediate() : curr_en(0), s_hash(0) {}

    // Constructor
    intermediate(size_t stringLength, size_t numMoves)
        : moves(numMoves, std::make_pair(0, 0)), // Initialize moves with numMoves number of (0,0) pairs
          curr_en(0), s_hash(0)                  // Initialize other members as needed
    {
    }
};

struct c_intermediate
{
    int curr_intermediate;
    int curr_en;
    int iter;
    int move_id;
    uint64_t s_hash;
};

void debug_print_pt(const short *temp_pt, int overrideNumberOfElements = -1)
{
    // prints a pairing table as if you would call print(pt) in Python
    if (temp_pt == nullptr)
    {
        std::cout << "[]" << std::endl;
        return;
    }

    int numberOfElements;
    if (overrideNumberOfElements != -1)
    {
        numberOfElements = overrideNumberOfElements + 1;
    }
    else
    {
        numberOfElements = temp_pt[0] + 1; // Including the first element itself
    }

    std::cout << "[";
    for (int i = 0; i < numberOfElements; ++i)
    {
        std::cout << temp_pt[i];
        if (i < numberOfElements - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

// Helper function to format the moves vector as a Python-like list
std::string format_moves(const std::vector<std::pair<int, int>> &moves)
{
    std::string formattedMoves = "[";
    for (size_t i = 0; i < moves.size(); ++i)
    {
        formattedMoves += "(" + std::to_string(moves[i].first) + ", " + std::to_string(moves[i].second) + ")";
        if (i < moves.size() - 1)
        {
            formattedMoves += ", ";
        }
    }
    formattedMoves += "]";
    return formattedMoves;
}

void update_pt(const auto &current_intermediate, short *pt)
{
    // Assuming vrna_db_from_ptable and formatMoves are defined elsewhere
    std::fill_n(pt + 1, pt[0], static_cast<short>(0));

    for (const auto &current_move : current_intermediate.moves)
    {
        if (current_move.first == 0)
        {
            break;
        }
        pt[current_move.first] = current_move.second;
        pt[current_move.second] = current_move.first;
    }
}

std::string intermediate_to_string(const auto &current_intermediate, short *pt)
{
    update_pt(current_intermediate, pt);
    std::string s = vrna_db_from_ptable(pt);
    return s;
}

void debug_print_intermediate(const auto &current_intermediate, short *pt)
{
    std::cout << intermediate_to_string(current_intermediate, pt) << "\t"
              << current_intermediate.s_hash << "\t"
              << format_moves(current_intermediate.moves) << "\t"
              << current_intermediate.curr_en << std::endl;
}

// Function to print the vector of intermediates
void debug_print_intermediates(const std::vector<intermediate> &intermediates, short *pt)
{
    // Print header
    std::cout << "Structure\tHash\tMoves\tCurrent Energy\n";

    // Print each intermediate's details
    for (const auto &inter : intermediates)
    {
        std::cout << intermediate_to_string(inter, pt) << "\t"
                  << inter.s_hash << "\t"
                  << format_moves(inter.moves) << "\t"
                  << inter.curr_en << std::endl;
    }
}

int predictNumBasePairs(int sequenceLength, double slope = 0.306620, double intercept = -3.622861, double correctionFactor = 1.2)
{
    double predictedBasePairs = slope * sequenceLength + intercept;
    predictedBasePairs *= correctionFactor;

    // std::ceil is used to round up to the nearest whole number
    return static_cast<int>(std::ceil(predictedBasePairs));
}

void loopidx_from_ptable(const short *pt, short *loop, short *stack)
{
    // computes loop tables without memory allocation
    int i, hx, l, nl;
    int length;

    length = pt[0];
    memset(loop + 1, 0, sizeof(short) * pt[0]);
    hx = l = nl = 0;

    for (i = 1; i <= length; i++)
    {
        if ((pt[i] != 0) && (i < pt[i]))
        {
            /* ( */
            nl++;
            l = nl;
            stack[hx++] = i;
        }

        loop[i] = l;

        if ((pt[i] != 0) && (i > pt[i]))
        {
            /* ) */
            --hx;
            if (hx > 0)
                l = loop[stack[hx - 1]]; /* index of enclosing loop   */
            else
                l = 0; /* external loop has index 0 */
        }
    }
    loop[0] = nl;
}

// Simplified version of fold function
std::string fold(const std::string &seq, float bp_threshold = 0.01, int bps = -1, int l = 10, bool debug = false)
{
    if (bps == -1)
    {
        bps = predictNumBasePairs(seq.length());
    }

    vrna_md_t md;
    vrna_fold_compound_t *fc;

    vrna_md_set_default(&md); // copy global settings
    fc = vrna_fold_compound(seq.c_str(), &md, VRNA_OPTION_EVAL_ONLY);

    // init structure string (for return string), pairing tables,
    // stack and loop table arrays for temporary storage
    std::string s(seq.length(), '.');
    short *temp_pt = vrna_ptable(s.c_str());
    short *temp_loop = vrna_ptable_copy(temp_pt);
    short *stack = vrna_ptable_copy(temp_pt);

    std::vector<std::pair<int, int>> moves;
    std::vector<uint64_t> hash_move_list{};

    if (bp_threshold == 0)
    {
        // convert the vrna moves (all possible base pairs) into C++
        // generate list of associated hash values as well.
        vrna_move_t *moves_c = vrna_neighbors(fc, temp_pt, VRNA_MOVESET_INSERTION);
        for (auto nb = moves_c; nb->pos_5 != 0 && nb->pos_3 != 0; nb++)
        {

            int i = nb->pos_5;
            int j = nb->pos_3;
            // printf("add %d %d\n", i, j);
            moves.emplace_back(i, j);
            uint64_t move_hash;
            move_hash = int_hash_64(-i);
            move_hash *= (-j);
            move_hash += (-j) * 13;
            hash_move_list.push_back(move_hash);
        }
        free(moves_c);
    }
    else
    {
        // use only moves corresponding to potential base pairs according to the ViennaRNA-computed
        // pairs with base pair probability above the selected threshold
        // https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/examples/c.html?highlight=base+pair+probabilities
        char *propensity = (char *)vrna_alloc(sizeof(char) * (seq.length() + 1));
        vrna_ep_t *ptr, *pair_probabilities = NULL;

        float en = vrna_pf_fold(seq.c_str(), propensity, &pair_probabilities);

        // printf("%s\n%s [ %6.2f ]\n", seq.c_str(), propensity, en);

        /* print all base pairs with probability above 50% */
        for (ptr = pair_probabilities; ptr->i != 0; ptr++)
        {
            if (ptr->p >= bp_threshold)
            {

                int i = ptr->i;
                int j = ptr->j;

                // printf("add %d %d\n", i, j);

                moves.emplace_back(i, j);
                uint64_t move_hash;
                move_hash = int_hash_64(-i);
                move_hash *= (-j);
                move_hash += (-j) * 13;
                hash_move_list.push_back(move_hash);
            }
        }
        // printf("p(%d, %d) = %g\n", ptr->i, ptr->j, ptr->p);
        /* cleanup memory */
        free(pair_probabilities);
        free(propensity);
    }

    // predict number of base pairs - this could need some improvement
    int num_bp = predictNumBasePairs(temp_pt[0]);
    std::vector<std::pair<int, int>> init_moves(num_bp, std::make_pair(0, 0));

    // init starting intermediates
    intermediate init_intermediate(temp_pt[0], num_bp);
    intermediate best_intermediate(temp_pt[0], num_bp);
    best_intermediate.curr_en = 99999;

    std::vector<intermediate> curr_intermediates;
    curr_intermediates.push_back(init_intermediate);

    // vector for offspring candidates per iteration, will be cleared in each iteration,
    // with reserved number of items for efficient memory allocations.
    std::vector<c_intermediate> candidates;
    const int max_elements = l * moves.size();
    candidates.reserve(max_elements);

    // this is the set which tracks duplicates per iteration of the offspring candidates.
    // it only stores the (future) hash values of the structure, estimated via integer
    // hashing. all this info is also inside the candidates vector, but it is more efficient
    // to store the hash values here separately.
    tsl::robin_set<int> next_hash_values = {};
    next_hash_values.reserve(max_elements);

    if (debug)
        debug_print_intermediates(curr_intermediates, temp_pt);

    for (int iter = 0; iter < init_moves.size(); ++iter)
    {

        candidates.clear();
        next_hash_values.clear();

        int curr_best_en = 99999;

        // std::unordered_set<int> next_hash_values = {};

        // iterate over retained moves
        for (int b = 0; b < curr_intermediates.size(); b++)
        {
            auto &current_intermediate = curr_intermediates[b];

            // both options work here
            update_pt(current_intermediate, temp_pt);

            // update temp_loop according to current temp_pt
            loopidx_from_ptable(temp_pt, temp_loop, stack);

            // iterate over the list of available moves
            for (int a = 0; a < moves.size(); a++)
            {
                auto curr_i = moves[a].first;
                auto curr_j = moves[a].second;
                // std::cout << "try " << curr_i << "/" << curr_j << std::endl;

                if (temp_pt[curr_i] == 0 and temp_pt[curr_j] == 0 and temp_loop[curr_i] == temp_loop[curr_j])
                {
                    // before any energy evaluation, check if this candidate isn't already known
                    // (duplicate check)
                    auto candidate_hash = hash_move_list[a] + current_intermediate.s_hash;
                    if (!next_hash_values.contains(candidate_hash))
                    {
                        next_hash_values.insert(candidate_hash);
                    }
                    else
                    {
                        continue;
                    }

                    int en = current_intermediate.curr_en + vrna_eval_move_pt(fc, temp_pt, curr_i, curr_j);
                    // en_candidate = fc.eval_move_pt(pt, i, j) + energy

                    // update structure hash value

                    // intermediate candidate_intermediate = current_intermediate;
                    // candidate_intermediate.moves[iter].first = curr_i;
                    // candidate_intermediate.moves[iter].second = curr_j;
                    // candidate_intermediate.curr_en = en;
                    // candidate_intermediate.s_hash = candidate_hash;
                    // next_intermediates.push_back(candidate_intermediate);

                    // c_intermediate test;

                    candidates.emplace_back(c_intermediate{b, en, iter, a, candidate_hash});

                    if (en < best_intermediate.curr_en)
                    {
                        // intermediate candidate_intermediate = current_intermediate;
                        // candidate_intermediate.moves[iter].first = curr_i;
                        // candidate_intermediate.moves[iter].second = curr_j;
                        // candidate_intermediate.curr_en = en;
                        // candidate_intermediate.s_hash = candidate_hash;
                        // best_intermediate = candidate_intermediate;
                        best_intermediate = current_intermediate;
                        best_intermediate.moves[iter].first = curr_i;
                        best_intermediate.moves[iter].second = curr_j;
                        best_intermediate.curr_en = en;
                        best_intermediate.s_hash = candidate_hash;
                    }
                    if (en < curr_best_en)
                    {
                        curr_best_en = en;
                    }

                    // std::cout << "add" << curr_i << "/" << curr_j << "/" << en << "/" << curr_hash << "\n";
                }
                else
                {
                    // std::cout << "invalid\n";
                }
            }

            // if (b == 2)
            // {
            //     break;
            // }
        }

        std::nth_element(candidates.begin(), candidates.begin() + l, candidates.end(),
                         [](const c_intermediate &a, const c_intermediate &b)
                         {
                             return a.curr_en < b.curr_en;
                         });
        // candidates.resize(l); // no need to resize, just copy the l best elements...

        std::vector<intermediate> next_intermediates;

        // int ll = 0;
        // for (const auto &inter : candidates)
        // {
        //     if (ll == l)
        //     {
        //         break;
        //     }
        //     intermediate candidate_intermediate = curr_intermediates[inter.curr_intermediate];
        //     candidate_intermediate.s_hash = inter.s_hash;
        //     candidate_intermediate.curr_en = inter.curr_en;
        //     const int i_move = moves[inter.move_id].first;
        //     const int j_move = moves[inter.move_id].second;
        //     candidate_intermediate.moves[iter].first = i_move;
        //     candidate_intermediate.moves[iter].second = j_move;
        //     next_intermediates.push_back(candidate_intermediate);
        //     ll++;
        // }

        for (int index = 0; index < l and index < candidates.size(); ++index)
        {
            const auto &inter = candidates[index];
            intermediate candidate_intermediate = curr_intermediates[inter.curr_intermediate];

            candidate_intermediate.s_hash = inter.s_hash;
            candidate_intermediate.curr_en = inter.curr_en;

            const int i_move = moves[inter.move_id].first;
            const int j_move = moves[inter.move_id].second;
            candidate_intermediate.moves[iter].first = i_move; // iter is from the overall iteration...
            candidate_intermediate.moves[iter].second = j_move;

            next_intermediates.push_back(candidate_intermediate);
        }

        // std::cout << "added:" << next_intermediates.size() << "/" << candidates.size() << "/" << l*moves.size() << "\n";

        if (debug)
        {
            debug_print_intermediates(next_intermediates, temp_pt);
            std::cout << "\n";
        }

        curr_intermediates = next_intermediates;

        // std::cout << "iteration " << iter << best_intermediate.curr_en << "/" << curr_best_en << "\n";

        // early termination
        // if (best_intermediate.curr_en < curr_best_en - 400)
        // {
        //     break;
        // }
    }

    if (debug)
    {
        std::cout << "best found intermediate: \n";
        debug_print_intermediate(best_intermediate, temp_pt);
    }

    // std::cout << "END\n " << std::endl;
    s = intermediate_to_string(best_intermediate, temp_pt);

    free(temp_pt);
    free(stack);
    free(temp_loop);

    return s;
}


int main(int argc, char *argv[])
{
    bool debug = false;
    // bool debug = true;

    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <sequence> <l> [bp_threshold]" << std::endl;
        return 1;
    }

    std::string seq = argv[1];
    int l = std::atoi(argv[2]);

    // Set default value for bp_threshold
    float bp_threshold = 0.01;

    // Check if bp_threshold is provided as an optional argument
    if (argc >= 4)
    {
        bp_threshold = std::atof(argv[3]);
    }

    // std::string folded = fold(seq, -1, l);
    std::cout << fold(seq, bp_threshold, -1, l, debug) << std::endl;

    return 0;
}