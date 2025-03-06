"""Implementation of continuity in evolution
Fontana, W., & Schuster, P. (1998). Continuity in evolution: on the nature of transitions. Science, 280(5368), 1451-1455.
"""

#### TODO (all done?)
# - Sequence population should be random at the beginning
# - Create way to evaluate RNA.tree_edit_distance(ss_algo, ss_RNAfold) < threshold.
# - Make sure LCA sequences are not dead.
# - Time should not be +1, but function of the fitness.
#####


# recent additions (17/1/2025)
#   _check_if_alive: checks if the predicted MFE from the vrna calculation matches
#                    the folding algorithm’s evaluated energy for that structure.
#
#   DeadReplicate Event: if the sequence is dead, it gets the label 'D', to effectively kill the replicate
#
#   class PathToRNAfold (Adaptive Walk to “Alive”)
#
#   new _setup to start with random sequences len(self.target) with RNA.random_string()


# new stuff (18/1/2025):
#   verbosity + print intervals
#
#   time rescaling: from an exerpimental point of view, if we would use time+=1 or time+=dt does
#   not change the final genetic outcome, only the timestamps on those events change.
#   Who replicates and which mutations occur are not influenced by the time increments,
#   but this changes the plotting!
#   
#   Notes about the current exponential time variant:
#       The distance‐vs‐time curve commonly “slows down” near the optimum.
#       Mathematically, that’s because there are many more replication events 
#       packed into the same or smaller “wall clock” period.
#       Meaning, if the population becomes collectively fitter, the expected 
#       waiting time per even gets smaller. As we get closer to the end, clock
#       increments happen more slowly from the simulation's point of view.
#       This is exactly what Gillespie-style simulations typically show. 


# Fix bug (22/1/2025)
# In ContinuousEvolBase._replace function, move target_pop.fold to location after the replicated sequence is validated (alive)

import math
import random
# import local lib
from .neutral import AdaptiveWalk
from pathlib import Path

import numpy as np

import RNA

# Resolution for replicate time
RepAcc = 999
NUC = set(list('ACGU'))
Mutation = {k: list(NUC-set(k)) for k in NUC}

def repNuc(nuc):
    """Replicate one nucleotide
    """
    r = random.randint(1, 1000)
    if r <= RepAcc:
        return nuc, False
    return random.choice(Mutation[nuc]), True


def repSeq(seq):
    """Replicate a sequence
    """
    seq, mutate = zip(*(x for x in map(repNuc, seq)))
    return ''.join(seq), any(mutate)


class Event:
    def __init__(self, label, time, pop_id):
        self.label = label
        self.time = time
        self.pop = pop_id

    def __iter__(self):
        return (self.__dict__[item] for item in ('time', 'pop', 'label'))

class Start(Event):
    def __init__(self, pop_id, init_seq, init_ss, folding):
        super().__init__('S', 0, pop_id)
        self.seq = init_seq
        self.ss = init_ss
        self.folding = folding

    def __iter__(self):
        return (self.__dict__[item] for item in ('time', 'pop', 'label', 'seq', 'ss', 'folding'))


class Break(Event):
    def __init__(self, time, pop_id, new_seq, new_ss, folding):
        super().__init__('B', time, pop_id)
        self.seq = new_seq
        self.ss = new_ss
        self.folding = folding

    def __iter__(self):
        return (self.__dict__[item] for item in ('time', 'pop', 'label', 'seq', 'ss', 'folding'))

class Replicate(Event):
    """Event class for replication
    """
    def __init__(self, time, pop_id, new_seq):
        super().__init__('R', time, pop_id)
        self.seq = new_seq

    def __iter__(self):
        return (self.__dict__[item] for item in ('time', 'pop', 'label', 'seq'))

class DeadReplicate(Event):
    """Event class for replication where replicate is dead
    """
    def __init__(self, time, pop_id, new_seq):
        super().__init__('D', time, pop_id)
        self.seq = new_seq

    def __iter__(self):
        return (self.__dict__[item] for item in ('time', 'pop', 'label', 'seq'))

class SingleEvol:
    """Class object for single sequence evolution
    """
    def __init__(self, target, folding, iden=None):
        self.target = target
        self.folding = folding
        self.id = iden
        self.targetLen = len(target)
        self.history = []
        self.current_seq = ""
        self.current_ss = ""

    def fold(self, seq, time, track=None, ss=None):
        """Fold a sequence at given time
        Note: calling this function implies a "break" in the evolution
            - A start of evolution if time is 0
            - The current population is replaced by a new one (given sequence) and folding algorithm might also change

        Args:
            seq: new sequence
            time: time in evolution
            ss: shortcut to avoid refolding if given (in the case where replicate is the same sequence)
        """
        if ss is None:
            ss = self.folding.fold(seq)
        if time == 0:
            event = Start(self.id, seq, ss, self.folding.name)
        else:
            event = Break(time, self.id, seq, ss, self.folding.name)
        # self.history.append(event)
        if track is not None:
            with open(track, 'a') as f:
                print(*event, sep='\t', file=f)
        self.current_seq = seq
        self.current_ss = ss
        return ss

    def replicate(self, time, track=None):
        new_seq, is_new = repSeq(self.current_seq)
        event = Replicate(time, self.id, new_seq)
        # self.history.append(event)
        if track is not None:
            with open(track, 'a') as f:
                print(*event, sep='\t', file=f)
        return new_seq, is_new

    def to_csv(self, path_to_file, sep='\t', include_replicate=False):
        with open(path_to_file, 'w') as f:
            for event in self.history:
                if include_replicate or (not isinstance(event, Replicate)):
                    print(*event, sep=sep, file=f)

def fitness(target, dbn):
    return 1/(0.01+RNA.bp_distance(dbn, target)/len(target))

class ContinuousEvolBase:
    """Basic class object to study population evolution
    """
    def __init__(self, nPop=1000, verbose=False, track=None, output="test"):
        """Initiation with given folding algorithm

        Args:
            target: target structure
            nPop: population size
            track: store entire evolution history in given file at each step. It's highly not recommended to turn on since it will slow down the calculation
        """
        # Initiate all variables
        self.verbose = verbose
        self.track = track
        self.target = ""
        self.init_seq = ""
        self.nPop = nPop
        self.population = []
        self.fitness_pop = []
        self.time = 0
        self.history = []
        self.output = Path(output)
        self.output.mkdir(exist_ok=True, parents=True)
        self.output_distance = Path(self.output / 'history.txt')

    def _force_dump(self):
        with self.output_distance.open('a') as f:
            print(*self.history, file=f)
        self.history = []

    def _update_history(self):
        pass

    def _setup(self):
        """Initial population and/or condition etc
        """
        pass

    def _stop_evol(self):
        """Return True if evolution should stop
        """
        return False

    def _chose_one_to_replicate(self):
        """Chose one individual from population for replication
        """
        return random.choices(self.population, weights=self.fitness_pop)[0]


    def _chose_one_to_be_replaced(self):
        """Chose one individual from population to remove
        """
        return random.choice(self.population)

    def _check_if_alive(self, folding, seq):
        fc = RNA.fold_compound(seq)
        _, mfe = fc.mfe()
        return mfe == fc.eval_structure(folding.fold(seq))


    def _replace(self, source_pop, target_pop):
        """Define how to replace one individual by another one (replicate)
        """
        new_seq, is_new = source_pop.replicate(self.time, track=self.track)
        ss = None if is_new else source_pop.current_ss # Avoid calculate ss if is old
        # Check if new replicate is still alive
        if is_new and not self._check_if_alive(source_pop.folding, new_seq):
            event = DeadReplicate(self.time, source_pop.id, new_seq)
            source_pop.history.append(event)
            return
        # fold can ONLY be called after checking the new sequence is still alive
        new_ss = target_pop.fold(new_seq, self.time, track=self.track, ss=ss)
        self.fitness_pop[target_pop.id] = fitness(self.target, new_ss) #  Recalculate fitness new_ss

        if self.verbose:
            self._print_replace_info(source_pop, target_pop)

        return 1

    def _print_replace_info(self, source_pop, target_pop):
        print("At t={}, {}th is replicated and {}th is replaced".format(self.time, source_pop.id, target_pop.id), end='\r')

    def _increment_time(self, parents_fitness):
        Lambda = sum(self.fitness_pop) # total rate Lambda = sum of all lambda_i (equal to fitness[i])
        U = random.random()
        dt = -math.log(U)/Lambda # dt from Exp(Lambda)

        dt *= 1000  # scale up times

        self.time += dt

        # previous approaches...
        # self.time += 1/parents_fitness
        # self.time += 1

        return

    def evolution(self, target):
        """Start evolution with giving initial sequence

        Args:
            target: target structure
            seq: initial sequence to evolve
        """
        # Init
        self.target = target
        # Configure all population
        self._setup()

        parent_fitness = 1

        # Start evolution
        # Stop evolution when strictly more that half of population reaches target
        # self.time = 1 # this is defined in the setup above?

        self.event_count = 0  # Discrete event counter

        while not self._stop_evol():
            # first increment the time 
            self._increment_time(parent_fitness)
            self.event_count += 1

            # Select one individual to replicate according to the fitness
            to_replicate = self._chose_one_to_replicate()
            parent_fitness = self.fitness_pop[to_replicate.id]
            # Select one in the old population to replace
            to_remove = self._chose_one_to_be_replaced()
            self._replace(to_replicate, to_remove)

            if self.verbose > 1 and self.event_count % self.verbosity_step == 0:
                print()  # end the progress bar line
                print(f"[Event {self.event_count}], time={self.time:.5f}:")
                for i, pop in enumerate(self.population[:3]):  # Show 3 samples
                    print(f"  {i}: Seq: {pop.current_seq}, SS: {pop.current_ss}, Dist: {RNA.bp_distance(self.target, pop.current_ss)}")
            else:
                # print(self.nbReached, end="\r")
                pass

    def dump_history(self, directory_to_store, sep='\t', include_replicate=False):
        """Dump evolution history into files.
        History of each individual in the population is stored in each single file
        """
        output = Path(directory_to_store)
        output.mkdir(exist_ok=True, parents=True)
        with Path(output / 'history.npy').open('wb') as f:
            np.save(f, self.history)
        for pop in self.population:
            pop.to_csv(Path(output/ f"history_pop_{pop.id}.csv"), sep=sep, include_replicate=include_replicate)

class PathToRNAfold(AdaptiveWalk):
    """Class object for a adaptive walk within sequence space till given folding algorithm predicts same as RNAfold
    """
    def init_state(self, trial_seq, **kwargs):
        super().init_state(trial_seq, **kwargs)

    def distance(self, cur_seq, **kwargs):
        if self.flag is None:
            return self.distance_bp(cur_seq, **kwargs)
        return self.distance_tree(cur_seq, **kwargs)

    def distance_tree(self, cur_seq, **kwargs):
        fc = RNA.fold_compound(cur_seq)
        ss = self.folding.fold(cur_seq)
        tree_ss = RNA.make_tree(RNA.db_to_tree_string(ss, self.flag))

        dist = []
        for x in fc.subopt(0):
            tree_mfe = RNA.make_tree(RNA.db_to_tree_string(x.structure, self.flag))
            dist.append(RNA.tree_edit_distance(tree_mfe, tree_ss))
        return min(dist)

    def distance_bp(self, cur_seq, **kwargs):
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

class ContinuousEvol(ContinuousEvolBase):
    """Class object to reproduce study continuity in evolution same folding algorithm for all population
    """

    def __init__(self, folding, nPop=1000, verbose=False, track=None, step=1, max_time=100000000, verbosity_step=100000000, output='test'):
        """Initiation with given folding algorithm

        Args:
            folding: folding algorithm
            target: target structure
            nPop: population size
        """
        super().__init__(nPop=nPop, verbose=verbose, track=track, output=output)
        self.folding = folding

        # Additional variables for stop condition
        self.step = step # save interval
        self.half_pop = self.nPop/2
        self.reached_pop = []
        self.distance = []
        self.nbReached = 0
        self.history = []
        #CAVH: time restriction to look_behind
        self.max_time = max_time
        #Initial bpdistance
        self.bpdistance = np.zeros((self.nPop+1, 2*self.nPop+1))
        #Initial treedistance
        self.treedist = np.zeros((self.nPop+1, 2*self.nPop+1))

        self.verbosity_step = verbosity_step


    def _setup(self):
        #initial list of random sequences. They are not guatanteed to be alive
        init_seq = []

        #Mutate the sequences such that they are alive
        while len(init_seq) < self.nPop:
            rand_seq = RNA.random_string(len(self.target),'ACGU')
            tmp = PathToRNAfold(self.folding)
            final = tmp(rand_seq, rounds=10, loose=True)
            if final.dist == 0:
                init_seq.append(final.seq)
        print("Initial alive Sequences generated")

        for ind,current_seq in enumerate(init_seq):
            s = SingleEvol(self.target, self.folding, ind) # s is for ind features
            s.fold(current_seq, 0, track=self.track)
            self.population.append(s)

            fit = fitness(self.target, s.current_ss)
            dist = RNA.bp_distance(self.target, s.current_ss)
            self.fitness_pop.append(fit)
            self.distance.append(dist)
            # self.history.append([0,dist])
            self.reached_pop.append(s.current_ss==self.target)
        self.nbReached = sum(self.reached_pop)
        self._update_history()

        if self.verbose > 1:
            print("Initial Population (5 members):")
            for i, pop in enumerate(self.population[:5]):
                print(f"  {i}: Seq: {pop.current_seq}, SS: {pop.current_ss}, Dist: {RNA.bp_distance(self.target, pop.current_ss)}")
            print ()

    def _update_history(self):
        self.history.append((self.time, sum(self.distance)/len(self.distance)))
        if len(self.history) >= 1000:
            self._force_dump()

    def _stop_evol(self):
        #return (self.nbReached > self.half_pop)
        return (self.nbReached > self.half_pop) or (self.time > self.max_time)

    def _replace(self, source_pop, target_pop):
        """Define how to replace one individual by another one (replicate)
        """
        status = super()._replace(source_pop, target_pop)
        # The replicate is dead
        if status is None:
            return

        # Update stop condition related variables
        new_ss = target_pop.current_ss
        # CAVH: comparison to RNAfold prediction on that seq. not target.
        self.distance[target_pop.id] = RNA.bp_distance(self.target, new_ss)
        # if self.time % self.step == 0: # QS: Intervals to summarise distance (mean)
        if self.event_count % self.step == 0:
            self._update_history()
        self.nbReached += (int(new_ss == self.target) - int(self.reached_pop[target_pop.id]))
        self.reached_pop[target_pop.id] = new_ss == self.target

        return 1

    def _print_replace_info(self, source_pop, target_pop):
        print("At t={}, {}th is replicated and {}th is replaced. Dist to target: {} # of reaches: {}".format(self.time, source_pop.id, target_pop.id, RNA.bp_distance(self.target, target_pop.current_ss), self.nbReached), end='\r')

    def parse(self):
        """Parse result of evolution
        """
        pass

    def evolution(self, target):
        super().evolution(target)

        if self.verbose > 1:
            print("\n\nFinal Population:")
            for i, pop in enumerate(self.population[:5]):  # Show 5 samples
                print(f"  {i}: Seq: {pop.current_seq}, SS: {pop.current_ss}, Dist: {RNA.bp_distance(self.target, pop.current_ss)}")

        # final save
        self._update_history()
        self._force_dump()

class CompetitiveEvol(ContinuousEvolBase):
    """Class object to compare two or more folding algorithms until one taking over in evolution
    """
    def __init__(self, foldings, weights, nPop=100, verbose=False, track=None, step=10, threshold=0.9, max_time=300000):
        assert len(foldings) == len(weights), "foldings and weights should have the same length"
        super().__init__(nPop=nPop, verbose=verbose, track=track)

        self.foldings = foldings
        self.weights = weights

        self.max_time = max_time
        self.step = step 
        self.indx = {folding.name: ind for ind, folding in enumerate(foldings)}
        self.counts = [nPop] * len(foldings)
        self.history = np.array([[0] + self.counts])
        self.nPop_threshold = threshold * nPop * len(foldings)

    def _setup(self):
        ind = 0
        for folding in self.foldings:
            for _ in range(self.nPop):
                s = SingleEvol(self.target, folding, ind)
                ss = s.fold(self.init_seq, 0, track=self.track)
                self.population.append(s)
                self.fitness_pop.append(fitness(self.target, ss) * self.weights[self.indx[folding.name]])
                ind += 1

    def _stop_evol(self):
        return (self.time > self.max_time) or any(x>=self.nPop_threshold for x in self.counts)

    def _replace(self, source_pop, target_pop):
        """Define how to replace one individual by another one (replicate)
        """
        before = target_pop.folding.name
        after = source_pop.folding.name
        new_seq, is_new = source_pop.replicate(self.time, track=self.track)
        target_pop.folding = source_pop.folding
        ss = None if is_new else source_pop.current_ss
        new_ss = target_pop.fold(new_seq, self.time, track=self.track, ss=ss)
        self.fitness_pop[target_pop.id] = fitness(self.target, new_ss) * self.weights[self.indx[after]]

        # Update stop condition related variables
        self.counts[self.indx[before]] -= 1
        self.counts[self.indx[after]] += 1

        # if self.time % self.step == 0:
        # if self.event_count % self.step == 0:
            # self.history = np.append(self.history, [[self.time] + self.counts], axis=0)

        if self.verbose:
            self._print_replace_info(source_pop, target_pop)

    def _print_replace_info(self, source_pop, target_pop):
        tmp = ', '.join(f'{fold.name}: {count}' for fold, count in zip(self.foldings, self.counts))
        print(f"At t={self.time}, {tmp}", end='\r')

    def evolution(self, target, init_seq):
        super().evolution(target, init_seq)
        self.history = np.append(self.history, [[self.time] + self.counts], axis=0)

