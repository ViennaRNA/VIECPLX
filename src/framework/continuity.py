"""Implementation of continuity in evolution
Fontana, W., & Schuster, P. (1998). Continuity in evolution: on the nature of transitions. Science, 280(5368), 1451-1455.
"""

import random
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
        self.history.append(event)
        if track is not None:
            with open(track, 'a') as f:
                print(*event, sep='\t', file=f)
        self.current_seq = seq
        self.current_ss = ss
        return ss

    def replicate(self, time, track=None):
        new_seq, is_new = repSeq(self.current_seq)
        event = Replicate(time, self.id, new_seq)
        self.history.append(event)
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
    def __init__(self, nPop=1000, verbose=False, track=None):
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
        self.history = np.array([[]])


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

    def _replace(self, source_pop, target_pop):
        """Define how to replace one individual by another one (replicate)
        """
        new_seq, is_new = source_pop.replicate(self.time, track=self.track)
        ss = None if is_new else source_pop.current_ss
        new_ss = target_pop.fold(new_seq, self.time, track=self.track, ss=ss)
        self.fitness_pop[target_pop.id] = fitness(self.target, new_ss)

        if self.verbose:
            self._print_replace_info(source_pop, target_pop)

    def _print_replace_info(self, source_pop, target_pop):
        print("At t={}, {}th is replicated and {}th is replaced".format(self.time, source_pop.id, target_pop.id), end='\r')

    def evolution(self, target, init_seq):
        """Start evolution with giving initial sequence

        Args:
            target: target structure
            seq: initial sequence to evolve
        """
        # Init
        self.target = target
        self.init_seq = init_seq
        self._setup()

        # Start evolution
        # Stop evolution when strictly more that half of population reaches target
        self.time = 1
        while not self._stop_evol():
            # Select one individual to replicate according to the fitness
            to_replicate = self._chose_one_to_replicate()
            # Select one in the old population to replace
            to_remove = self._chose_one_to_be_replaced()
            self._replace(to_replicate, to_remove)

            self.time += 1

    def dump_history(self, directory_to_store, sep='\t', include_replicate=False):
        """Dump evolution history into files.
        Hisotry of each individual in the population is stored in each single file
        """
        output = Path(directory_to_store)
        output.mkdir(exist_ok=True, parents=True)
        with Path(output / 'history.npy').open('wb') as f:
            np.save(f, self.history)
        for pop in self.population:
            pop.to_csv(Path(output/ f"history_pop_{pop.id}.csv"), sep=sep, include_replicate=include_replicate)


class ContinuousEvol(ContinuousEvolBase):
    """Class object to reproduce study continuity in evolution same folding algorithm for all population
    """

    def __init__(self, folding, nPop=1000, verbose=False, track=None, step=1, max_time=100000000):
        """Initiation with given folding algorithm

        Args:
            folding: folding algorithm
            target: target structure
            nPop: population size
        """
        super().__init__(nPop=nPop, verbose=verbose, track=track)
        self.folding = folding

        # Additional variables for stop condition
        self.step = step
        self.half_pop = self.nPop/2
        self.reached_pop = []
        self.distance = []
        self.nbReached = 0
        #CAVH: time restriction to look_behind
        self.max_time = max_time

    def _setup(self):
        ss = self.folding.fold(self.init_seq)

        for ind in range(self.nPop):
            s = SingleEvol(self.target, self.folding, ind)
            s.fold(self.init_seq, 0, track=self.track, ss=ss)
            self.population.append(s)

        fit = fitness(self.target, ss)
        dist = RNA.bp_distance(self.target, ss)
        self.fitness_pop = [fit] * self.nPop
        self.distance = np.array([dist] * self.nPop)
        self.history = np.array([[0, dist]])
        self.reached_pop = [ss==self.target] * self.nPop
        self.nbReached = sum(self.reached_pop)

    def _stop_evol(self):
        #return (self.nbReached > self.half_pop)
        return (self.nbReached > self.half_pop) or (self.time > self.max_time)

    def _replace(self, source_pop, target_pop):
        """Define how to replace one individual by another one (replicate)
        """
        super()._replace(source_pop, target_pop)

        # Update stop condition related variables
        new_ss = target_pop.current_ss
        self.distance[target_pop.id] = RNA.bp_distance(self.target, new_ss)
        if self.time % self.step == 0:
            self.history = np.append(self.history, [[self.time, np.mean(self.distance)]], axis=0)
        self.nbReached += (int(new_ss == self.target) - int(self.reached_pop[target_pop.id]))
        self.reached_pop[target_pop.id] = new_ss == self.target

    def _print_replace_info(self, source_pop, target_pop):
        print("At t={}, {}th is replicated and {}th is replaced. Dist to target: {} # of reaches: {}".format(self.time, source_pop.id, target_pop.id, RNA.bp_distance(self.target, target_pop.current_ss), self.nbReached), end='\r')

    def parse(self):
        """Parse result of evolution
        """
        pass

    def evolution(self, target, init_seq):
        super().evolution(target, init_seq)
        self.history = np.append(self.history, [[self.time, np.mean(self.distance)]], axis=0)

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

        if self.time % self.step == 0:
            self.history = np.append(self.history, [[self.time] + self.counts], axis=0)

        if self.verbose:
            self._print_replace_info(source_pop, target_pop)

    def _print_replace_info(self, source_pop, target_pop):
        tmp = ', '.join(f'{fold.name}: {count}' for fold, count in zip(self.foldings, self.counts))
        print(f"At t={self.time}, {tmp}", end='\r')

    def evolution(self, target, init_seq):
        super().evolution(target, init_seq)
        self.history = np.append(self.history, [[self.time] + self.counts], axis=0)

