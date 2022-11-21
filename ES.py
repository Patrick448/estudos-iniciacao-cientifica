import datetime
import math
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


def alg1plus1ES(expr="", validator=None, sigma=1.0, c=0.817, n=10, iter=100, seed=0):
    np.random.seed(seed)
    p = []
    A = []
    t = 0
    ind = {}
    obj_func_hist = []

    for v in parser.parse(expr).variables():
        ind[v] = 0

    p.append(ind)
    phi = evaluate(expr, validator, **ind)
    calculated_error = None

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

        p.append(mutated_ind)
        obj_func_hist.append(new_phi)
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

        # if t % (20*n) == 0 and t != 0:
        #     calculated_error = abs(evaluate(expr, validator, **ind) - evaluate(expr, validator, **p[len(p) - 20*n]))
        #     if abs(evaluate(expr, validator, **ind) - evaluate(expr, validator, **p[len(p) - 20 * n])) <= error:
        #         break

    return p, p[-2:-1], obj_func_hist


def populational_isotropic_ES(expr="", validator=None, dimension_gen_interval=(0, 0), sigma_var=0.5, iter=100, seed=0,
                              num_parents=0, num_offspring=0):
    np.random.seed(seed)
    population = []
    t = 0

    for i in range(num_parents):
        ind = {'dim': {}, 'sigma': None}
        for v in parser.parse(expr).variables():
            ind['dim'][v] = np.random.uniform(dimension_gen_interval[0], dimension_gen_interval[1], 1)[0]
        ind['sigma'] = np.random.uniform(0, 1, 1)[0]
        phi = evaluate(expr, validator, **ind['dim'])
        ind['eval'] = phi
        population.append(ind)

    # todo: mudar criterio para chamadas da função objetivo talvez

    for i in range(iter):
        t += 1
        offspring = []

        for parent in population:
            for j in range(math.ceil(num_offspring / num_parents)):
                mutated_ind = {'dim': {}, 'sigma': None}
                new_sigma = parent['sigma'] * math.exp(np.random.normal() * sigma_var ** 2)  # muta o sigma
                mutated_ind['sigma'] = new_sigma
                for key in parent['dim']:  # muta cada uma das dimensões
                    mutated_ind['dim'][key] = parent['dim'][key] + new_sigma * np.random.normal()

                mutated_ind['eval'] = evaluate(expr, validator, **mutated_ind['dim'])  # faz avaliação
                offspring.append(mutated_ind)

        offspring_and_parents = population + offspring
        best = sorted(offspring_and_parents, key=lambda i: i['eval'])[0:num_parents]
        population = best

    return population


def populational_non_isotropic_ES(expr="", validator=None, dimension_gen_interval=(0, 0), sigma_var=0.5, iter=100, seed=0,
                              num_parents=0, num_offspring=0):
    np.random.seed(seed)
    population = []
    t = 0

    for i in range(num_parents):
        ind = {'dim': {}, 'sigma': {}}
        for v in parser.parse(expr).variables():
            ind['dim'][v] = np.random.uniform(dimension_gen_interval[0], dimension_gen_interval[1], 1)[0]
            ind['sigma'][v] = np.random.uniform(0, 1, 1)[0]
        phi = evaluate(expr, validator, **ind['dim'])
        ind['eval'] = phi
        population.append(ind)

    # todo: mudar criterio para chamadas da função objetivo talvez

    for i in range(iter):
        t += 1
        offspring = []

        for parent in population:
            for j in range(math.ceil(num_offspring / num_parents)):
                mutated_ind = {'dim': {}, 'sigma': {}}
                for key in parent['dim']:  # muta cada uma das dimensões e seus sigmas
                    mutated_ind['sigma'][key] = parent['sigma'][key] * math.exp(np.random.normal() * sigma_var ** 2)
                    mutated_ind['dim'][key] = parent['dim'][key] + mutated_ind['sigma'][key] * np.random.normal()

                mutated_ind['eval'] = evaluate(expr, validator, **mutated_ind['dim'])  # faz avaliação
                offspring.append(mutated_ind)

        offspring_and_parents = population + offspring
        best = sorted(offspring_and_parents, key=lambda i: i['eval'])[0:num_parents]
        population = best

    return population


# alg1plus1ES(1, 10, 0.817, 10, 100, int(time.time()))


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


class Individual:
    sigmas = {}
    dimensions = {}
    age = 0
    evaluation = 0

    def __init__(self):
        pass

    def __gt__(self, other):
        return self.get_evaluation() > other.get_evaluation()

    def __lt__(self, other):
        return self.get_evaluation() < other.get_evaluation()

    def __str__(self):
        return str(self.to_dict())

    def to_dict(self):
        return {'dim': self.dimensions, 'sigma': self.sigmas, 'eval': self.evaluation}

    def get_evaluation(self):
        return self.evaluation

    def set_evaluation(self, evaluation):
        self.evaluation = evaluation

    def add_dimension(self, name, value, sigma):
        self.dimensions[name] = value
        self.sigmas[name] = sigma

    def set_info(self, dim_dict, sigmas_dict, eval):
        self.dimensions = dim_dict
        self.sigmas = sigmas_dict
        self.evaluation = eval

    def set_dimensions(self, dim_dict):
        self.dimensions = dim_dict

    def set_sigmas(self, sigmas_dict):
        self.sigmas = sigmas_dict

    def get_dimensions(self):
        return self.dimensions

    def get_sigmas(self):
        return self.sigmas


class ESAlgorithm:
    evaluation_function = None
    interval_validator = None
    population = []
    population_history = []
    success_history = []

    def __init__(self):
        self.evaluation_function = None
        self.interval_validator = None
        self.population = []
        self.population_history = []
        self.success_history = []
        pass

    def set_evaluation_function(self, func):
        self.evaluation_function = func

    def set_interval_validator(self, func):
        self.interval_validator = func

    def simple_mutation(self, individual, sigma):
        dims = individual.get_dimensions()
        new_ind = Individual()

        for key in dims:
            dim_value = dims[key] + sigma * np.random.normal()
            new_ind.add_dimension(key, dim_value, sigma)

        new_ind.set_evaluation(self.evaluation_function(new_ind.get_dimensions()))
        print(f"Anterior: {individual}")
        print(f"Novo: {new_ind}")
        return new_ind

    def run_one_plus_one(self, first_ind={}, sigma=1.0, c=0.817, n=10, iterations=100, seed=0):
        np.random.seed(seed)
        t = 0
        current_sigma = sigma
        sigmas = {}
        for key in first_ind:
            sigmas[key] = current_sigma

        ind = Individual()
        ind.set_info(first_ind.copy(), sigmas, self.evaluation_function(first_ind))

        self.population.append(ind)
        calculated_error = None

        for i in range(iterations):
            t += 1
            mutated_ind = self.simple_mutation(ind, current_sigma)

            if mutated_ind < ind:
                ind = mutated_ind
                self.success_history.append(1)
            else:
                self.success_history.append(0)

            self.population.append(mutated_ind)

            if t % n == 0:
                # get successes and failures from at most 10n entries in A
                window = self.success_history[-10 * n:]
                success = sum(window)
                failure = abs(success - len(window))
                ps = success / (success + failure)
                if ps < 1 / 5:
                    current_sigma = current_sigma * c
                elif ps > 1 / 5:
                    current_sigma = current_sigma / c

        return self.population, ind
