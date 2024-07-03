ranks = [2, 3, 5, 10, 20, 30, 50]
degrees = [2, 3, 5, 10, 20, 30, 50]
algorithms = ["motive", "crystalline", "CSA"]

fp = sys.stdout
Fq = GF(5)
A.<T> = Fq[]
heads = ["rank", "[K:Fq]"] + algorithms
print("\t".join(heads))
fp.flush()
for r in ranks:
    for d in degrees:
        row = [str(r), str(d)]
        K.<z> = Fq.extension(d)
        total = walltime()
        timings = { algorithm: 0 for algorithm in algorithms }
        repeat = 0
        while walltime(total) < 10:
            coeffs = [z] + [K.random_element() for _ in range(r-2)] + [1]
            phi = DrinfeldModule(A, coeffs)
            for algorithm in algorithms:
                t = walltime()
                _ = phi.frobenius_charpoly(algorithm=algorithm)
                timings[algorithm] += walltime(t)
            repeat += 1
        for algorithm in algorithms:
            row.append("%.6fs" % (timings[algorithm]/repeat))
        row.append("(repeat=%s)" % repeat)
        print("\t".join(row))
        fp.flush()
