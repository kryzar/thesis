algorithms = ["motive", "crystalline", "CSA"]
timeout = 10

#r = 10
#degrees = [2, 3, 5, 10, 20, 30, 50, 100]
d = 15
ranks = [3, 5, 10, 20, 30, 50, 100]

yscale=0.5
xscale=1.2
pad = 0.05
opacity = 0.8

fp = sys.stdout

Fq = GF(5)
A.<T> = Fq[]
lines = { algorithm: [ ] for algorithm in algorithms }
#for d in degrees:
for r in ranks:
    K.<z> = Fq.extension(d)
    total = walltime()
    timings = { algorithm: 0 for algorithm in algorithms }
    repeat = 0
    while walltime(total) < timeout:
        coeffs = [z] + [K.random_element() for _ in range(r-2)] + [1]
        phi = DrinfeldModule(A, coeffs)
        for algorithm in algorithms:
            t = walltime()
            _ = phi.frobenius_charpoly(algorithm=algorithm)
            timings[algorithm] += walltime(t)
        repeat += 1
    #row = [str(d)]
    row = [str(r)]
    for algorithm in algorithms:
        row.append("%.6fs" % (timings[algorithm]/repeat))
        #x = xscale * float(log(d))
        x = xscale * float(log(r))
        y = yscale * log(timings[algorithm] / repeat)
        lines[algorithm].append((x,y))
    row.append("(repeat=%s)" % repeat)
    print("% " + "\t".join(row))
    fp.flush()

print("\\begin{tikzpicture}")
for algorithm in algorithms:
    print("\\begin{scope}[%s, opacity=%.1f, transparency group]" % (algorithm, opacity))
    path = [ "(%.3f,%.3f)" % (x, y) for x, y in lines[algorithm] ]
    print("  \draw " + "--".join(path) + ";")
    for x, y in lines[algorithm]:
        print("  \\fill (%.3f,%.3f) rectangle (%.3f,%.3f);" % (x-pad, y-pad, x+pad, y+pad))
    print("\\end{scope}")
print("\\end{tikzpicture}")
