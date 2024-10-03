
# A Relation between Natural Selection and Computational Complexity

Source code and script to reproduce the figures presented in the paper

## Project description

The space of possible structures that a given RNA sequence can form is large ($2.7^n$). If the RNA (biological system) would have to explore all possible structures to find the correct functional structure it would take astronomical time rather than being compatible with timeframes of biological processes.

While computational algorithms for counting or finding the optimal structure are NP-hard(for certain definitions...), biological systems evolved to solve the problem much faster - for a subset of sequences. Algorithms on a secondary structure level can find optimal structures within $O(n^3)$, which is most likely (?) still significantly slower than what biological systems implement through evolution. Can we find algorithms with lower complexity that still solve the folding problem - possibly for subsets of RNA sequences?



## Code and script

- `src/` contains the source code of presented
  - folding algorithms (`foldingAlg`), including look behind folding, basic co-fold, best helix co-fold (named Best helix cofold with RNAfold in second step in the source code), folding rule, and beam search.
  - simulation experiments (`framework`), including structure density surface (sds), shape frequency, upper bound neutral path, and evolution.
- `script/` contains the script to run each sumulation experiment of folding algorithms. The result is stored in `result/`
- `doc/` contains jupyter notebook that plots the figures presented in the paper from results stored in `result/`
