from pathlib import Path
from datetime import datetime
from time import time
from statistics import median
import logging


################
# BASIC CONFIG #
################


logging.basicConfig(format="[%(asctime)s] %(message)s", level=logging.INFO)
set_random_seed(4808)


#####################
# GLOBAL PARAMETERS #
#####################


# EXTENSION_DEGREES = [2, 3, 10, 50, 100, 500, 1000]
# TAU_DEGREES       = [2, 3, 10, 50, 100, 500, 1000] 
# RANKS             = [2, 3, 10, 50, 100, 500, 1000]

# DEFAULT_EXTENSION_DEGREE   = 15
# DEFAULT_RANK               = 10
# DEFAULT_ISOGENY_TAU_DEGREE = 10

EXTENSION_DEGREES = [2, 3]
TAU_DEGREES       = [2, 3] 
RANKS             = [2, 3]

DEFAULT_EXTENSION_DEGREE   = 15
DEFAULT_RANK               = 10
DEFAULT_ISOGENY_TAU_DEGREE = 10

NUMBER_SAMPLES = 2


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
    """
    Rapidly find an endomorphism of `phi` whose tau-degree is lesser
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
        endomorphism *= u * v + w
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
            isogeny = ore_pol * isogeny
            assert not phi.is_isomorphic(psi)
            # logging.info(isogeny.degree())
        except:
            pass
    # logging.info('Excited while')
    isogeny = Hom(phi, psi)(isogeny)
    # logging.info('Created the isogeny object')
    return isogeny


##########################
# BENCHMARKING FUNCTIONS #
##########################


def get_samples(f, phi, n, r, d, param, is_isogeny):
    iso_or_endo      = 'isogeny' if is_isogeny else 'endomorphism'
    norm_or_charpoly = 'norm'    if is_isogeny else 'charpoly'
    # Get samples
    samples = []
    for sample_number in range(NUMBER_SAMPLES):
        logging.info(f'(n, r, d) = ({n}, {r}, {d})')
        logging.info(f'Param: {param}')
        logging.info(f'Sample: {sample_number+1}/{NUMBER_SAMPLES}')
        logging.info(f'Starting {iso_or_endo} computation...')
        morphism = find_isogeny(phi, n) if is_isogeny else find_endomorphism(phi, n)
        logging.info(f'{iso_or_endo} computed.')
        callable = morphism.norm if is_isogeny else morphism.charpoly
        logging.info(f'Starting {iso_or_endo} {norm_or_charpoly} computation...')
        computation_time = time_callable(callable)
        samples.append(computation_time)
        logging.info(f'Isogeny norm computed ({computation_time}s).')
    # Write to a file
    match param:
        case 'n':
           abscissa = n
        case 'd':
           abscissa = d
        case 'r':
           abscissa = r
    f.write(f'{abscissa} {median(samples)}\n')

def bench_tau_degree(filename, is_isogeny):
    r = DEFAULT_RANK
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    drinfeld_modules = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        for n in TAU_DEGREES:
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, 'n', is_isogeny)

def bench_extension_degree(filename, is_isogeny):
    r = DEFAULT_RANK
    n = DEFAULT_ISOGENY_TAU_DEGREE
    # Get computation times
    with open(filename, 'w') as f:
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
        for r in RANKS:
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, 'r', is_isogeny)


if __name__ == '__main__':

    Fq = GF(5)
    A.<T> = Fq[]

    workdir = Path.home() / Path('workshop/benchmarks/tests')
    date = datetime.now().strftime('%H:%M:%S')

    for bench_function in [bench_tau_degree, bench_extension_degree, bench_rank]:
        for is_isogeny in [True, False]:
            norm_or_charpoly = 'norm' if is_isogeny else 'charpoly'
            file = Path(f'{date}' \
                        f'.{norm_or_charpoly}' \
                        f'.{bench_function.__name__}' \
                        f'.txt')
            bench_function(workdir / file, is_isogeny)

