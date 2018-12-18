"""
pyart.lazydict
==============

A dictionary-like class supporting lazy loading of specified keys.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    LazyLoadDict

"""

try:
    # Python 3
    from collections.abc import MutableMapping
except ImportError:
    # Python 2.7, will be removed in next release after Py-ART Impressionism.
    from collections import MutableMapping
import itertools


class LazyLoadDict(MutableMapping):
    """
    A dictionary-like class supporting lazy loading of specified keys.

    Keys which are lazy loaded are specified using the set_lazy method.
    The callable object which produces the specified key is provided as the
    second argument to this method.  This object gets called when the value
    of the key is loaded. After this initial call the results is cached
    in the traditional dictionary which is used for supplemental access to
    this key.

    Testing for keys in this dictionary using the "key in d" syntax will
    result in the loading of a lazy key, use "key in d.keys()" to prevent
    this evaluation.

    The comparison methods, __cmp__, __ge__, __gt__, __le__, __lt__, __ne__,
    nor the view methods, viewitems, viewkeys, viewvalues, are implemented.
    Neither is the the fromkeys method.

    Parameters
    ----------
    dic : dict
        Dictionary containing key, value pairs which will be stored and
        evaluated traditionally. This dictionary referenced not copied into
        the LazyLoadDictionary and hence changed to this dictionary may change
        the original.  If this behavior is not desired copy dic in the
        initalization.

    Examples
    --------
    >>> d = LazyLoadDict({'key1': 'value1', 'key2': 'value2'})
    >>> d.keys()
    ['key2', 'key1']
    >>> lazy_func = lambda : 999
    >>> d.set_lazy('lazykey1', lazy_func)
    >>> d.keys()
    ['key2', 'key1', 'lazykey1']
    >>> d['lazykey1']
    999

    """

    def __init__(self, dic):
        """ initalize. """
        self._dic = dic
        self._lazyload = {}

    # abstract methods
    def __setitem__(self, key, value):
        """ Set a key which will not be stored and evaluated traditionally. """
        self._dic[key] = value
        if key in self._lazyload:
            del self._lazyload[key]

    def __getitem__(self, key):
        """ Get the value of a key, evaluating a lazy key if needed. """
        if key in self._lazyload:
            value = self._lazyload[key]()
            self._dic[key] = value
            del self._lazyload[key]
        return self._dic[key]

    def __delitem__(self, key):
        """ Remove a lazy or traditional key from the dictionary. """
        if key in self._lazyload:
            del self._lazyload[key]
        else:
            del self._dic[key]

    def __iter__(self):
        """ Iterate over all lazy and traditional keys. """
        return itertools.chain(self._dic.copy(), self._lazyload.copy())

    def __len__(self):
        """ Return the number of traditional and lazy keys. """
        return len(self._dic) + len(self._lazyload)

    # additional class to mimic dict behavior
    def __str__(self):
        """ Return a string representation of the object. """
        if len(self._dic) == 0 or len(self._lazyload) == 0:
            seperator = ''
        else:
            seperator = ', '
        lazy_reprs = [(repr(k), repr(v)) for k, v in self._lazyload.items()]
        lazy_strs = ['%s: LazyLoad(%s)' % r for r in lazy_reprs]
        lazy_str = ", ".join(lazy_strs) + '}'
        return str(self._dic)[:-1] + seperator + lazy_str

    def has_key(self, key):
        """ True if dictionary has key, else False. """
        return key in self

    def copy(self):
        """
        Return a copy of the dictionary.

        Lazy keys are not evaluated in the original or copied dictionary.
        """
        dic = self.__class__(self._dic.copy())
        # load all lazy keys into the copy
        for key, value_callable in self._lazyload.items():
            dic.set_lazy(key, value_callable)
        return dic

    # lazy dictionary specific methods
    def set_lazy(self, key, value_callable):
        """ Set a lazy key to load from a callable object. """
        if key in self._dic:
            del self._dic[key]
        self._lazyload[key] = value_callable
