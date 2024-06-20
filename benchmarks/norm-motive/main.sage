from datetime import datetime
from time import time
from statistics import median
import logging

logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.INFO)

EXTENSION_DEGREES   = [2, 3, 5, 10, 20, 30, 50, 100]
ISOGENY_TAU_DEGREES = [2, 3, 5, 10, 20, 30, 50, 100] 
RANKS               = [2, 3, 5, 10, 20, 30, 50, 100]

DEFAULT_EXTENSION_DEGREE = 15
DEFAULT_RANK             = 10
DEFAULT_TAU_DEGREE       = 10

NUMBER_SAMPLES = 20


def time_callable(callable):
    """
    Return the time require to call `callable`, which is typically
    a function, e.g. `phi.frobenius_charpoly`.
    """
    start = time()
    callable()
    end = time()
    return end - start


def find_isogenous(phi):
    tau_K = phi.frobenius_endomorphism().ore_polynomial()
    while True:
        (u, v) = (A.random_element(), A.random_element())
        isogeny = phi(u).right_gcd(tau_K - phi(v))
        try:
            psi = phi.velu(isogeny)
            if not phi.is_isomorphic(psi):
                return psi
        except:
            pass


def non_zero_isogeny(phi, psi, n):
    while True:
        morphism = Hom(phi, psi).random_element(n)
        if morphism.is_isogeny():
             return morphism


def bench_tau_degree(filename):
    r = DEFAULT_RANK
    d = DEFAULT_EXTENSION_DEGREE
    K.<z> = Fq.extension(d)
    category = DrinfeldModule(A, [z, 1]).category()
    # Get computation times
    with open(filename, 'w') as f:
        for n in ISOGENY_TAU_DEGREES:
            phi = category.random_object(r)
            psi = find_isogenous(phi)
            samples = []
            for sample_number in range(NUMBER_SAMPLES):
                logging.warning(f'n={n}, sample {sample_number+1}/{NUMBER_SAMPLES}')
                isogeny = non_zero_isogeny(phi, psi, n)
                samples.append(time_callable(isogeny.norm))
            # Process samples
            median_comptime = median(samples)
            logging.info(f'median comptime: {median_comptime}')
            # Write them to a file
            f.write(f'{n}    {median_comptime}\n')


if __name__ == '__main__':

    Fq = GF(5)
    A.<T> = Fq[]
    filename = '/users/ahugoune/workshop/benchmarks/norm_charpoly-tau_degree.' \
               + datetime.now().strftime('%Y-%m-%dT%H:%M:%S%Z') \
               + '.txt'
    bench_tau_degree(filename)
