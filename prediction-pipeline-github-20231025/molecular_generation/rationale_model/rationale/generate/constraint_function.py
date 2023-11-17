from functools import partial

def tddft_func(x):
    return min(x) >= 500

def sa_func(x):
    return min(x) <= 3.5

def tddft_sa_func(x):
    return float(x[0]) >= float(500) and float(x[1]) <= float(3.5)

def normal_func(x):
    return min(x) >= 0.5

def tddft_ensemble_old(x):
    return min(x) >= 200

def constraint_func_old(x):
    return tddft_ensemble(x)

def tddft_ensemble(cutoff,x):
    return min(x) >= cutoff

def tddft_ensemble_func(cutoff):
    return partial(tddft_ensemble,cutoff)

