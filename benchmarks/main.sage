from pathlib import Path
from datetime import datetime
from time import time
from statistics import median
import logging

logging.basicConfig(format="[%(asctime)s:%(funcName)s()] %(message)s",
			        level=logging.INFO)

set_random_seed(4808)

EXTENSION_DEGREES = [2, 5, 10, 50, 100, 500, 1000]
TAU_DEGREES       = [2, 5, 10, 50, 100, 500, 1000] 
RANKS             = [2, 5, 10, 50, 100, 500, 1000]

DEFAULT_EXTENSION_DEGREE   = 15
DEFAULT_RANK               = 10
DEFAULT_ISOGENY_TAU_DEGREE = 10

NUMBER_SAMPLES = 10


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
    end_ = End(phi)
    if n <= 4:
        endomorphism = end_.random_element(n)
    else:
        endomorphism = end_(1)
        basis_4 = end_.Fq_basis(4)
        basis_5 = end_.Fq_basis(5)
        basis_6 = end_.Fq_basis(6)
        Fq = phi.function_ring().base()
        print('Hey')
        while endomorphism.ore_polynomial().degree() < n:
            print(endomorphism.ore_polynomial().degree())
            u = sum((Fq.random_element() * x for x in basis_4))
            v = sum((Fq.random_element() * x for x in basis_5))
            w = sum((Fq.random_element() * x for x in basis_6))
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


def get_samples(f, phi, n, r, d, norm_or_charpoly):
    samples = []
    # Get samples
    for sample_number in range(NUMBER_SAMPLES):
        logging.info(f'(n, r, d) = ({n}, {r}, {d})')
        logging.info('sample {sample_number+1} out of {NUMBER_SAMPLES}')
        match norm_or_charpoly:
            case 'norm':
                logging.info('Starting isogeny computation...')
                isogeny = find_isogeny(phi, n)
                logging.info('Isogeny computed.')
                logging.info('Starting isogeny norm computation...')
                samples.append(time_callable(isogeny.norm))
                logging.info('Isogeny norm computed.')
            case 'charpoly':
                logging.info('Starting charpoly computation...')
                endomorphism = find_charpoly(phi, n)
                logging.info('Isogeny computed.')
                logging.info('Starting isogeny norm computation...')
                samples.append(time_callable(endomorphism.charpoly))
                logging.info('Isogeny norm computed.')


def bench_tau_degree(filename, norm_or_charpoly):
    r = DEFAULT_RANK
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    drinfeld_modules = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        for n in TAU_DEGREES:
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, norm_or_charpoly)


def bench_extension_degree(filename, norm_or_charpoly):
    r = DEFAULT_RANK
    n = DEFAULT_ISOGENY_TAU_DEGREE
    # Get computation times
    with open(filename, 'w') as f:
        for d in EXTENSION_DEGREES:
            K.<z> = Fq.extension(d)
            drinfeld_modules = DrinfeldModule(A, [z, 1]).category()
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, norm_or_charpoly)


def bench_rank(filename):
    n = DEFAULT_ISOGENY_TAU_DEGREE
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    drinfeld_modules = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        for r in RANKS:
            phi = drinfeld_modules.random_object(r)
            get_samples(f, phi, n, r, d, norm_or_charpoly)


if __name__ == '__main__':

    # Math
    Fq = GF(5)
    A.<T> = Fq[]
    # Boilerplate
    workdir = Path.home() / Path('workshop/benchmarks')
    date = datetime.now().strftime('%m-%d_%H:%M')
    for norm_or_charpoly in ['norm', 'charpoly']:
        bench_tau_degree(      norm_or_charpoly,
                               f'{workdir}/{date}.{norm_or_charpoly}.tau_degree.txt')
        bench_extension_degree(norm_or_charpoly,
                               f'{workdir}/{date}.{norm_or_charpoly}.extension_degree.txt')
        bench_rank(            norm_or_charpoly,
                               f'{workdir}/{date}.{norm_or_charpoly}.rank.txt')
