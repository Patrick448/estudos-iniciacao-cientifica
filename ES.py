import datetime
import random
import time
import numpy as np
from py_expression_eval import Parser

parser = Parser()


def f(x):
    return x * x + 6 * x - 4


def evaluate(expr, in_allowed_inteval, **kwargs):
    if in_allowed_inteval(**kwargs):
        return parser.parse(expr).evaluate(kwargs)

    return float('inf')


def in_allowed_interval_g(**kwargs):
    if -1.5 <= kwargs['x'] <= 4 and -3.0 <= kwargs['y'] <= 4:
        return True

    return False


def alg1plus1ES(expr, validator, sigma=1.0, c=0.817, n=10, iter=100, seed=0):
    np.random.seed(seed)
    p = []
    A = []
    t = 0
    ind = {}

    for v in parser.parse(expr).variables():
        ind[v] = 0

    p.append(ind)
    phi = evaluate(expr, validator, **ind)

    # todo: ver outros criterios de parada
    for i in range(iter):
        t += 1
        mutated_ind = {}

        # mutate individual
        for key in ind:
            mutated_ind[key] = ind[key] + sigma * np.random.normal()

        # calculate new phi based on the mutation
        new_phi = evaluate(expr, validator, **mutated_ind)

        if new_phi < phi:
            ind = mutated_ind
            phi = new_phi
            A.append(1)
            pass
        else:
            A.append(0)
            pass

        p.append(ind)
        if t % n == 0:
            # get successes and failures from at most 10n entries in A
            window = A[-10 * n:]
            success = sum(window)
            failure = abs(success - len(window))
            ps = success / (success + failure)
            if ps < 1 / 5:
                sigma = sigma * c
                pass
            elif ps > 1 / 5:
                sigma = sigma / c
                pass

    return p, p[-2:-1]


# alg1plus1ES(1, 10, 0.817, 10, 100, int(time.time()))


class Individual:
    def __init__(self, age, x, eval):
        self.age = age
        self.eval = eval
        self.x = x


def a(m, l):
    p = []
    for i in range(m):
        ind = Individual(1, 0, f(0))
        p.append(ind)

    while True:
        q = []

        for i in range(l):
            # randomly select p parents in p[]
            q_ind = Individual(0, 0, f(0))  # generate q from variation
            q.append(q_ind)

        # p = select the best m individuals from q + p
        # increase age by one of p's individuals and update psi.v
