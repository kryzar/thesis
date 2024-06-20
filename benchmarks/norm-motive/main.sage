from time import time
import logging

logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.INFO)

EXTENSION_DEGREES   = [2, 3, 5, 10, 20, 30, 50, 100]
ISOGENY_TAU_DEGREES = [2, 3, 5, 10, 20, 30, 50, 100] 
RANKS               = [2, 3, 5, 10, 20, 30, 50, 100]

DEFAULT_EXTENSION_DEGREE = 15
DEFAULT_RANK             = 10
DEFAULT_TAU_DEGREE       = 10

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


def dict_to_file(filename, dict_):
    with open(filename, 'w') as f:
        for key, value in dict_.items():
            f.write(f'{key}    {value}\n')


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
    computation_times = dict()
    for n in ISOGENY_TAU_DEGREES:
        phi = category.random_object(r)
        psi = find_isogenous(phi)
        computation_time = 0
        for sample in range(NUMBER_SAMPLES):
            logging.warning(f'n={n}, sample {sample+1}/{NUMBER_SAMPLES}')
            isogeny = non_zero_isogeny(phi, psi, n)
            computation_time += time_callable(isogeny.norm)
        computation_time /= NUMBER_SAMPLES
        logging.info(f'computation time: {computation_time}')
        computation_times[n] = computation_time
    # Write them to a file
    dict_to_file(filename, computation_times)


if __name__ == '__main__':

    Fq = GF(5)
    A.<T> = Fq[]
    bench_tau_degree('/users/ahugoune/workshop/benchmarks/norm_charpoly-tau_degree.txt')
