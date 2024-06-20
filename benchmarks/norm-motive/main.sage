from time import time

EXTENSION_DEGREES   = [2, 3, 5, 10, 20, 30, 50, 100]
ISOGENY_TAU_DEGREES = [2, 3, 5, 10, 20, 30, 50, 100] 
RANKS               = [2, 3, 5, 10, 20, 30, 50, 100]

DEFAULT_EXTENSION_DEGREE = 15
DEFAULT_RANK             = 10
DEFAULT_TAU_DEGREE       = 10

NUMBER_SAMPLES = 10

Fq = GF(5)
A.<T> = Fq[]


def time_callable(callable):
    """
    Return the time require to call `callable`, which is typically
    a function, e.g. `phi.frobenius_charpoly`.
    """
    start = time()
    _ = callable()
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
    with open(filename, 'r') as f:
        for key, value in dict_.items():
            f.write(f'{key}    {value}')


def bench_tau_degree(filename):
    """
    Bench.
    """
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
        for _ in NUMBER_SAMPLES:
            isogeny = Hom(phi, psi).random_element(n)
            computation_time += time_callable(isogeny.norm)
        computation_times[n] = computation_time / NUMBER_SAMPLES
    # Write them to a file
    dict_to_file(filename, computation_times)


if __name__ == '__main__':

    bench_tau_degree(filename)
