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


def time_callable(callable):
    """
    Return the time require to call `callable`, which is typically
    a function, e.g. `phi.frobenius_charpoly`.
    """
    start = time()
    callable()
    end = time()
    return end - start


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


def bench_tau_degree(filename):
    r = DEFAULT_RANK
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    category = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        for n in TAU_DEGREES:
            phi = category.random_object(r)
            samples = []
            for sample_number in range(NUMBER_SAMPLES):
                logging.info(f'n={n}, sample {sample_number+1}/{NUMBER_SAMPLES}')
                isogeny = find_isogeny(phi, n)
                logging.info('Isogeny computed')
                samples.append(time_callable(isogeny.norm))
                logging.info('Isogeny norm computed')
            # Process samples
            median_comptime = median(samples)
            logging.info(f'median comptime: {median_comptime}')
            # Write them to a file
            f.write(f'{n} {median_comptime}\n')


def bench_extension_degree(filename):
    r = DEFAULT_RANK
    n = DEFAULT_ISOGENY_TAU_DEGREE
    # Get computation times
    with open(filename, 'w') as f:
        for d in EXTENSION_DEGREES:
            K.<z> = Fq.extension(d)
            category = DrinfeldModule(A, [z, 1]).category()
            phi = category.random_object(r)
            samples = []
            for sample_number in range(NUMBER_SAMPLES):
                logging.info(f'n={d}, sample {sample_number+1}/{NUMBER_SAMPLES}')
                isogeny = find_isogeny(phi, n)
                logging.info('Isogeny computed')
                samples.append(time_callable(isogeny.norm))
                logging.info('Isogeny norm computed')
            # Process samples
            median_comptime = median(samples)
            logging.info(f'median comptime: {median_comptime}')
            # Write them to a file
            f.write(f'{d} {median_comptime}\n')


def bench_rank(filename):
    n = DEFAULT_ISOGENY_TAU_DEGREE
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    category = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        for r in RANKS:
            phi = category.random_object(r)
            samples = []
            for sample_number in range(NUMBER_SAMPLES):
                logging.info(f'n={r}, sample {sample_number+1}/{NUMBER_SAMPLES}')
                isogeny = find_isogeny(phi, n)
                logging.info('Isogeny computed')
                samples.append(time_callable(isogeny.norm))
                logging.info('Isogeny norm computed')
            # Process samples
            median_comptime = median(samples)
            logging.info(f'median comptime: {median_comptime}')
            # Write them to a file
            f.write(f'{r} {median_comptime}\n')


if __name__ == '__main__':

    Fq = GF(5)
    A.<T> = Fq[]
    workdir = '/users/ahugoune/workshop/benchmarks'
    date = datetime.now().strftime('%Y-%m-%dT%H:%M:%S%Z')
    bench_tau_degree(      f'{workdir}/{date}.norm_charpoly-tau_degree.txt')
    bench_extension_degree(f'{workdir}/{date}.norm_charpoly-extension_degree.txt')
    bench_rank(            f'{workdir}/{date}.norm_charpoly-rank.txt')
