import RNA


__all__ = ['Folding', 'ViennaFold', 'Vienna']


class Folding:
    """Basic class of folding algorithm
    """

    def __init__(self, name, algorithm=None):
        """Init folding object with given name and folding algorithm

        Args:
            name (str): name of folding algorithm
            algorithm: function returning structure in dbn for given sequence
        """
        self.name = name
        self._fold = algorithm

    def fold(self, w):
        """Return structure for given sequence
        """
        if self._fold is not None:
            return self._fold(w)
        else:
            raise TypeError('Folding algorithm is undefined')

# There are two ways to create a Folding object that will can be used in the anlysis
# First, we can initiate the basic object Folding with folding algorithm in the form of function
# For the use of multiprocessing later, one needs to avoid of using lambda


def _RNAfold(w):
    return RNA.fold(w)[0]


ViennaFold      = Folding('RNAfold', _RNAfold)
# Second, one can also create a new Folding class inheriting from Folding
# In such way, one will need to overwrite member function fold


class Vienna(Folding):
    def fold(self, w):
        """Overwrite for ViennaRNA
        """
        return RNA.fold(w)[0]
