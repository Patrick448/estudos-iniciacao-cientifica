import random

import numpy as np

np.random.seed(1)
print(np.random.normal())


def f(x):
    return x*x +6*x -4


def alg1plus1ES():
    p = [123.0]
    phi = f(p[0])
    A = []
    t = 0
    n = 10
    sigma = 1.0
    x = p[0]
    c = 0.817

    # todo: ver outros criterios de parada
    for i in range(1000):
        t += 1
        x2 = x + sigma*np.random.normal()
        phi2 = f(x2)
        if phi2 < phi:
            x = x2
            phi = phi2
            A.append(1)
            pass
        else:
            A.append(0)
            pass

        p.append(x)
        if t % n == 0:
            #get successes and failures from at most 10n entries in A
            window = A[-10*n:]
            success = sum(window)
            failure = abs(success - len(window))
            ps = success/(success+failure)
            if ps < 1 / 5:
                sigma = sigma*c
                pass
            elif ps > 1 / 5:
                sigma = sigma/c
                pass
    print("\n".join(str(x) for x in p))


alg1plus1ES()