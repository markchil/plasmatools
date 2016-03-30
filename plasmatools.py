from __future__ import division

import MDSplus
import scipy
import warnings

def shot_length(shot=None, tree='magnetics', threshold=1.2e5):
    """Gets the termination time of the plasma current.
    
    Based on the IDL program shot_length.pro.
    
    Valid calling patterns::
        
        # Use shot 1140729030, default to "magnetics" tree:
        l = shot_length(1140729030)
        
        # Same as above but with explicit keyword:
        l = shot_length(shot=1140729030)
        
        # Use shot 1140729030, but get Ip from the "cmod" tree:
        l = shot_length(1140729030, 'cmod')
        
        # Same as above but with explicit keywords:
        l = shot_length(shot=1140729030, tree='cmod')
        
        # Use a tree which has already been opened (saves overhead):
        t = MDSplus.Tree(1140729030, 'magnetics')
        l = shot_length(tree=t)
        
    
    Parameters
    ----------
    shot : int, optional
        The shot number to use. If absent, a :py:class:`MDSplus.Tree` instance
        must be given in the `tree` keyword.
    tree : str or :py:class:`MDSplus.Tree` instance, optional
        The tree to fetch the plasma current data from. If given as a string,
        `shot` must also be provided. The specified tree will be opened, and
        "\\" + `tree` + "::Ip" will be retrieved. If given as an
        :py:class:`MDSplus.Tree`, there must be a node named "Ip" at the top
        level. "magnetics" and "cmod" are good options.
    threshold : float, optional
        The plasma current threshold (in amperes) for plasma termination.
        Default is 1.2e5 (120 kA).
    """
    try:
        if isinstance(tree, MDSplus.Tree):
            d = str(tree.getDefault())
            try:
                n = tree.getNode(d[:-3] + "Ip")
            except:
                n = tree.getNode(d + "::Ip")
        else:
            tr = MDSplus.Tree(tree, shot)
            n = tr.getNode("\\" + tree + "::Ip")
        Ip = scipy.absolute(scipy.asarray(n.data(), dtype=float))
        t = scipy.asarray(n.dim_of().data(), dtype=float)
    except:
        # If the signals cannot be obtained, return 0.0 but print a warning:
        warnings.warn("Could not load plasma current data!", RuntimeWarning)
        return 0.0
    
    # Mask the times which are >= 0:
    mask = t >= 0.0
    Ip = Ip[mask]
    t = t[mask]
    # Form the Ip - thresh difference array:
    delta_prod_mask = scipy.concatenate((
        ((Ip[:-1] - threshold) * (Ip[1:] - threshold)) <= 0.0,
        Ip[-1:] <= threshold
    ))
    # Form the Ip difference arrays:
    delta_1_mask = scipy.concatenate((
        (Ip[1:] - Ip[:-1]) <= 0.0,
        Ip[-1:] <= threshold
    ))
    delta_8_mask = scipy.concatenate((
        (Ip[8:] - Ip[:-8]) <= 0.0,
        Ip[-8:] <= threshold
    ))
    # Find the first point where the conditions are satisfied:
    cond_met, = scipy.nonzero(delta_prod_mask & delta_1_mask & delta_8_mask)
    if len(cond_met) == 0:
        return 0.0
    else:
        return t[cond_met[0]]
