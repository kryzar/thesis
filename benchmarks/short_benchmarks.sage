from datetime import datetime
from hashlib import md5
from itertools import product
from logging import FileHandler, Formatter, StreamHandler, getLogger, INFO
from multiprocessing import Pool
from pathlib import Path
from statistics import mean, median, stdev
from time import time


###########
# LOGGING #
###########


DATE = datetime.now().strftime('%d-%m-%Y_%H:%M:%S')
WORKDIR = Path.home() / Path('workshop/benchmarks/')

logger = getLogger()
logger.setLevel(INFO)
file_handler = FileHandler(WORKDIR / Path(f'{DATE}.log'))
console_handler = StreamHandler()
handlers = [file_handler, console_handler]
for handler in handlers:
    handler.setLevel(INFO)
    handler.setFormatter(Formatter('[%(asctime)s] %(message)s'))
    logger.addHandler(handler)

def hash_id(var):
    hasher = md5()
    hasher.update(var.encode())
    return hasher.hexdigest()[0:8]

#####################
# GLOBAL PARAMETERS #
#####################


set_random_seed(4808)

RANGES =   list(range(1, 10)) \
         + list(range(10, 20, 2)) \
         + list(range(20, 50, 5)) \
         + list(range(50, 100, 10))

EXTENSION_DEGREES = RANGES
TAU_DEGREES       = RANGES
RANKS             = RANGES

DEFAULT_EXTENSION_DEGREE   = 15
DEFAULT_RANK               = 10
DEFAULT_ISOGENY_TAU_DEGREE = 10

NUMBER_SAMPLES = 20


############################
# MEASURE COMPUTATION TIME #
############################


def time_callable(callable):
    """
    Return the time require to call `callable`, which is typically
    a function, e.g. `phi.frobenius_charpoly`.
    """
    start = time()
    callable()
    end = time()
    return end - start


################################
# FIND ISOGENY OR ENDOMORPHISM #
################################


def find_endomorphism(phi, n):
    """ Rapidly find an endomorphism of `phi` whose tau-degree is lesser
    than `n`.
    """
    end_ = End(phi)
    Fq = phi.function_ring().base()
    # Compute bases of endomorphisms
    min_deg = min(phi.rank(),
                  phi.frobenius_endomorphism().ore_polynomial().degree())
    basis_u = end_.Fq_basis(min_deg)
    basis_v = end_.Fq_basis(min_deg + 1)
    basis_w = end_.Fq_basis(min_deg + 2)
    # Iteratively take random combinations of those endomorphisms:
    endomorphism = end_(1)
    while endomorphism.ore_polynomial().degree() < n:
        u = sum((Fq.random_element() * x for x in basis_u))
        v = sum((Fq.random_element() * x for x in basis_v))
        w = sum((Fq.random_element() * x for x in basis_w))
        u_v_w = u * v + w
        if not u_v_w.is_zero():
            endomorphism *= u_v_w
    return endomorphism


def find_isogeny(phi, n):
    """
    Being given a Drinfeld module `phi` and an integer `n`, find an
    isogeny from `phi` to another Drinfeld module that is not
    isomorphic to `phi`.
    """

    tau_K = phi.frobenius_endomorphism().ore_polynomial()
    isogeny = phi.ore_polring()(1)
    psi = phi
    while isogeny.degree() < n:
        try:
            (u, v) = (A.random_element(), A.random_element())
            ore_pol = psi(u).right_gcd(tau_K - psi(v))
            psi = psi.velu(ore_pol)
            if not ore_pol.is_zero():
                isogeny = ore_pol * isogeny
            assert not phi.is_isomorphic(psi)
        except:
            pass
    isogeny = Hom(phi, psi)(isogeny)
    return isogeny


##########################
# BENCHMARKING FUNCTIONS #
##########################


def get_samples(f, phi, n, r, d, param, is_isogeny):
    iso_or_endo      = 'isogeny' if is_isogeny else 'endomorphism'
    norm_or_charpoly = 'norm'    if is_isogeny else 'charpoly'
    # Get samples
    samples = []
    hash_ = hash_id(repr(locals()))
    for sample_number in range(NUMBER_SAMPLES):
        logger.info(f'[{hash_}] (n, r, d) = ({n}, {r}, {d})')
        logger.info(f'[{hash_}] Param: {param}')
        logger.info(f'[{hash_}] Sample: {sample_number+1}/{NUMBER_SAMPLES}')
        logger.info(f'[{hash_}] Starting {iso_or_endo} computation...')
        morphism = find_isogeny(phi, n) if is_isogeny else find_endomorphism(phi, n)
        logger.info(f'{iso_or_endo} computed.')
        callable = morphism.norm if is_isogeny else morphism.charpoly
        logger.info(f'[{hash_}] Starting {iso_or_endo} {norm_or_charpoly} computation...')
        computation_time = time_callable(callable)
        samples.append(computation_time)
        logger.info(f'[{hash_}] {iso_or_endo} {norm_or_charpoly} computed ({computation_time}s).')
    # Write to a file
    match param:
        case 'n':
           abscissa = n
        case 'd':
           abscissa = d
        case 'r':
           abscissa = r
    data = f'{abscissa} '        \
           f'{min(samples)} '    \
           f'{max(samples)} '    \
           f'{mean(samples)} '   \
           f'{median(samples)} ' \
           f'{stdev(samples)}'
    logger.info(f'[{hash_}] data [(n, d, r) = ({n}, {d}, {r}), param: {param}]: {data}')
    f.write(data + '\n')

def bench_tau_degree(filename, is_isogeny):
    r = DEFAULT_RANK
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    drinfeld_modules = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        f.write('abscissa min max mean median stdev\n')
        for n in TAU_DEGREES:
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, 'n', is_isogeny)

def bench_extension_degree(filename, is_isogeny):
    r = DEFAULT_RANK
    n = DEFAULT_ISOGENY_TAU_DEGREE
    # Get computation times
    with open(filename, 'w') as f:
        f.write('abscissa min max mean median stdev\n')
        for d in EXTENSION_DEGREES:
            K.<z> = Fq.extension(d)
            drinfeld_modules = DrinfeldModule(A, [z, 1]).category()
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, 'd', is_isogeny)

def bench_rank(filename, is_isogeny):
    n = DEFAULT_ISOGENY_TAU_DEGREE
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    drinfeld_modules = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        f.write('abscissa min max mean median stdev\n')
        for r in RANKS:
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, 'r', is_isogeny)


###################
# PARALLELIZATION #
###################


def call_bench_function(args):
    (bench_function, is_isogeny) = args
    norm_or_charpoly = 'norm' if is_isogeny else 'charpoly'
    file = Path(f'{DATE}.{norm_or_charpoly}.{bench_function.__name__}.txt')
    bench_function(WORKDIR / file, is_isogeny)


if __name__ == '__main__':

    Fq = GF(5)
    A.<T> = Fq[]
    bench_functions = [bench_tau_degree, bench_extension_degree, bench_rank]
    # Parallelization: one core for each bench function/input:
    with Pool() as pool:
        pool.map(call_bench_function, product(bench_functions, [True, False]))
