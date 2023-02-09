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
            p.append(mutated_ind)
            obj_func_hist.append(new_phi)
        else:
            A.append(0)

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
                new_sigma = parent['sigma'] * math.exp(np.random.normal(0, sigma_var ** 2))  # muta o sigma
                mutated_ind['sigma'] = new_sigma
                for key in parent['dim']:  # muta cada uma das dimensões
                    mutated_ind['dim'][key] = parent['dim'][key] + np.random.normal(0, new_sigma)

                mutated_ind['eval'] = evaluate(expr, validator, **mutated_ind['dim'])  # faz avaliação
                offspring.append(mutated_ind)

        offspring_and_parents = population + offspring
        best = sorted(offspring_and_parents, key=lambda i: i['eval'])[0:num_parents]
        population = best

    return population


def populational_non_isotropic_ES(expr="", validator=None, dimension_gen_interval=(0, 0), sigma_var=0.5, iter=100,
                                  seed=0,
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
                    mutated_ind['sigma'][key] = parent['sigma'][key] * math.exp(
                        np.random.normal(0, sigma_var ** 2)) * math.exp(np.random.normal(0, sigma_var ** 2))
                    mutated_ind['dim'][key] = parent['dim'][key] + np.random.normal(0, mutated_ind['sigma'][key])

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
    evaluation_expression = None
    evaluation_function = None
    interval_validator = None
    population = []
    population_history = []
    success_history = []
    variable_bounds = {}
    variable_bounds_test = []
    using_global_bounds = False
    global_variable_bounds = {}
    dimension_mapping = []
    num_dimensions = 0
    parser_expression_evaluator = None
    evaluation_counter = 0
    known_minimum = None
    execution_history = {'checkpoints':[], 'values':[]}

    def __init__(self):
        self.error_stop_criterion = 0
        self.using_stop_criterion = False
        self.evaluation_expression = None
        self.evaluation_function = None
        self.interval_validator = None
        self.population = []
        self.population_history = []
        self.success_history = []
        self.variable_bounds = {}
        self.using_global_bounds = False
        self.global_variable_bounds = {'upper':
                                           {'value': float('inf'), 'closed': True},
                                       'lower':
                                           {'value': float('-inf'), 'closed': True}
                                       }
        self.variable_bounds_test = []
        self.num_dimensions = 0
        self.execution_history = {'checkpoints': [], 'values': []}
        self.known_minimum = None

    def set_evaluation_function(self, func):
        self.evaluation_function = func

    def reset_execution_history(self):
        self.execution_history = {'checkpoints': [], 'values': []}

    def add_checkpoint(self, checkpoint, value):
        self.execution_history['checkpoints'].append(checkpoint)
        self.execution_history['values'].append(value)

    def default_expression_evaluator(self, ind):
        pass

    def set_evaluation_expression(self, expression):
        self.evaluation_expression = expression
        self.parser_expression_evaluator = parser.parse(self.evaluation_expression).evaluate

    def set_interval_validator(self, expression):
        self.interval_validator = expression

    def set_dimension_mapping(self, mapping):
        self.dimension_mapping = mapping
        self.num_dimensions = len(mapping)

    def set_num_dimensions(self, num):
        self.num_dimensions = num
        self.variable_bounds_test = [None] * self.num_dimensions

    def set_error_stop_criterion(self, error):
        self.using_stop_criterion = True
        self.error_stop_criterion = error

    def set_known_minimum(self, value):
        self.known_minimum = value

    def abs_error(self, val):
        if self.known_minimum:
            return abs(val - self.known_minimum)

        return 0

    def stop(self, val):
        if self.using_stop_criterion and self.abs_error(val) <= self.error_stop_criterion:
            return True

        return False

    def calculate_checkpoint(self, k, max_iter):
        return self.num_dimensions**(k/5 - 3)*max_iter

    def get_execution_history(self):
        return self.execution_history

    def get_evaluations_count(self):
        return self.evaluation_counter

    def set_global_variable_bounds(self, upper_bound=float('inf'), upper_bound_closed=True, lower_bound=float('-inf'),
                                   lower_bound_closed=True):
        self.using_global_bounds = True
        self.global_variable_bounds = {'upper':
                                           {'value': upper_bound, 'closed': upper_bound_closed},
                                       'lower':
                                           {'value': lower_bound, 'closed': lower_bound_closed}
                                       }

    def set_variable_bounds(self, variable, upper_bound=float('inf'), upper_bound_closed=True,
                            lower_bound=float('-inf'), lower_bound_closed=True):
        self.using_global_bounds = False
        self.variable_bounds[variable] = {'upper':
                                              {'value': upper_bound, 'closed': upper_bound_closed},
                                          'lower':
                                              {'value': lower_bound, 'closed': lower_bound_closed}}

    def set_variable_bounds_test(self, variable, upper_bound=float('inf'), upper_bound_closed=True,
                                 lower_bound=float('-inf'), lower_bound_closed=True):
        self.variable_bounds_test[variable] = {'upper':
                                                   {'value': upper_bound, 'closed': upper_bound_closed},
                                               'lower':
                                                   {'value': lower_bound, 'closed': lower_bound_closed}}

    def evaluate(self, expr, **kwargs):
        return parser.parse(expr).evaluate(kwargs)

    def evaluate_test(self, ind):
        self.evaluation_counter += 1

        if self.evaluation_expression is not None:
            dim_dict = {}
            for i, key in enumerate(self.dimension_mapping):
                dim_dict[key] = ind['dim'][i]
            return self.parser_expression_evaluator(dim_dict)
        elif self.evaluation_function is not None:
            return self.evaluation_function(ind['dim'])

    def get_bound_value(self, dimension, bound="upper"):

        if self.using_global_bounds:
            return self.global_variable_bounds[bound]['value']
        else:
            return self.variable_bounds[dimension][bound]['value']

    def get_bound_value_test(self, dimension, bound="upper"):

        if self.using_global_bounds:
            return self.global_variable_bounds[bound]['value']
        else:
            return self.variable_bounds_test[dimension][bound]['value']


    def validate(self, individual):

        # todo: check if global bounds are enabled, if true validate using it. if false, validade like before
        validated_individual = individual
        for key in individual['dim']:
            if self.variable_bounds[key]['upper']['closed'] and individual['dim'][key] > \
                    self.variable_bounds[key]['upper']['value']:
                validated_individual['dim'][key] = self.variable_bounds[key]['upper']['value']
            elif not self.variable_bounds[key]['upper']['closed'] and individual['dim'][key] >= \
                    self.variable_bounds[key]['upper']['value']:
                validated_individual['dim'][key] = self.variable_bounds[key]['upper']['value'] - 0.0000001
            elif self.variable_bounds[key]['lower']['closed'] and individual['dim'][key] < \
                    self.variable_bounds[key]['lower']['value']:
                validated_individual['dim'][key] = self.variable_bounds[key]['lower']['value']
            elif not self.variable_bounds[key]['lower']['closed'] and individual['dim'][key] <= \
                    self.variable_bounds[key]['lower']['value']:
                validated_individual['dim'][key] = self.variable_bounds[key]['lower']['value'] + 0.0000001

        return validated_individual

    def validate_test(self, individual):

        # todo: check if global bounds are enabled, if true validate using it. if false, validade like before
        validated_individual = individual

        upper_bound = self.global_variable_bounds['upper']['value']
        lower_bound = self.global_variable_bounds['lower']['value']
        upper_bound_closed = self.global_variable_bounds['upper']['closed']
        lower_bound_closed = self.global_variable_bounds['lower']['closed']
        for i, value in enumerate(individual['dim']):

            if not self.using_global_bounds:
                upper_bound = self.variable_bounds_test[i]['upper']['value']
                lower_bound = self.variable_bounds_test[i]['lower']['value']
                upper_bound_closed = self.variable_bounds_test[i]['upper']['closed']
                lower_bound_closed = self.variable_bounds_test[i]['lower']['closed']

            if upper_bound_closed and value > upper_bound:
                validated_individual['dim'][i] = upper_bound
            elif not upper_bound_closed and value >= upper_bound:
                validated_individual['dim'][i] = upper_bound - 0.0000001
            elif lower_bound_closed and value < lower_bound:
                validated_individual['dim'][i] = lower_bound
            elif not lower_bound_closed and value <= lower_bound:
                validated_individual['dim'][i] = lower_bound + 0.0000001

        return validated_individual

    def one_plus_one_ES(self, sigma=1.0, c=0.817, n=10, iter=100, seed=0):
        np.random.seed(seed)
        p = []
        A = []
        t = 0
        ind = {'dim': {}}
        obj_func_hist = []

        for v in parser.parse(self.evaluation_expression).variables():
            ind['dim'][v] = 0

        ind = self.validate(ind)
        ind['eval'] = self.evaluate(self.evaluation_expression, **ind['dim'])
        p.append(ind)

        for i in range(iter):
            t += 1
            mutated_ind = {'dim': {}}

            # mutate individual
            for key in ind['dim']:
                mutated_ind['dim'][key] = ind['dim'][key] + sigma * np.random.normal()

            mutated_ind = self.validate(mutated_ind)
            # calculate new phi based on the mutation
            mutated_ind['eval'] = self.evaluate(self.evaluation_expression, **mutated_ind['dim'])

            if mutated_ind['eval'] < ind['eval']:
                ind = mutated_ind
                A.append(1)
                p.append(mutated_ind)
            else:
                A.append(0)

            if t % n == 0:
                # get successes and failures from at most 10n entries in A
                window = A[-10 * n:]
                success = sum(window)
                failure = abs(success - len(window))
                ps = success / (success + failure)
                if ps < 1 / 5:
                    sigma = sigma * c
                elif ps > 1 / 5:
                    sigma = sigma / c



        return p

    def populational_isotropic_ES(self, dimension_gen_interval=(0, 0), sigma_var=0.5, iter=100, seed=0,
                                  num_parents=0, num_offspring=0):
        np.random.seed(seed)
        population = []
        t = 0

        for i in range(num_parents):
            individual = {'dim': {}, 'sigma': None}
            for v in parser.parse(self.evaluation_expression).variables():
                individual['dim'][v] = np.random.uniform(dimension_gen_interval[0], dimension_gen_interval[1], 1)[0]
            individual['sigma'] = np.random.uniform(0, 1, 1)[0]
            individual = self.validate(individual)
            evaluation = self.evaluate(self.evaluation_expression, **individual['dim'])
            individual['eval'] = evaluation
            population.append(individual)

        # todo: mudar criterio para chamadas da função objetivo talvez

        for i in range(iter):
            t += 1
            offspring = []

            for parent in population:
                for j in range(math.ceil(num_offspring / num_parents)):
                    mutated_ind = {'dim': {}, 'sigma': None}
                    mutated_ind['sigma'] = parent['sigma'] * math.exp(
                        np.random.normal(0, sigma_var ** 2))  # muta o sigma
                    for key in parent['dim']:  # muta cada uma das dimensões e seus sigmas
                        mutated_ind['dim'][key] = parent['dim'][key] + np.random.normal(0, mutated_ind['sigma'])
                    mutated_ind = self.validate(mutated_ind)
                    mutated_ind['eval'] = self.evaluate(self.evaluation_expression,
                                                        **mutated_ind['dim'])  # faz avaliação
                    offspring.append(mutated_ind)

            offspring_and_parents = population + offspring
            best = sorted(offspring_and_parents, key=lambda i: i['eval'])
            population = best[0:num_parents]

        return population

    def populational_non_isotropic_ES(self, dimension_gen_interval=(0, 0), sigma_var=0.5, iter=100, seed=0,
                                      num_parents=0, num_offspring=0):
        np.random.seed(seed)
        population = []
        t = 0

        for i in range(num_parents):
            individual = {'dim': {}, 'sigma': {}}
            for v in parser.parse(self.evaluation_expression).variables():
                individual['dim'][v] = np.random.uniform(dimension_gen_interval[0], dimension_gen_interval[1], 1)[0]
                individual['sigma'][v] = np.random.uniform(0, 1, 1)[0]
            individual = self.validate(individual)
            evaluation = self.evaluate(self.evaluation_expression, **individual['dim'])
            individual['eval'] = evaluation
            population.append(individual)

        # todo: mudar criterio para chamadas da função objetivo talvez

        for i in range(iter):
            t += 1
            offspring = []

            for parent in population:
                for j in range(math.ceil(num_offspring / num_parents)):
                    mutated_ind = {'dim': {}, 'sigma': {}}
                    for key in parent['dim']:  # muta cada uma das dimensões e seus sigmas
                        mutated_ind['sigma'][key] = parent['sigma'][key] * math.exp(
                            np.random.normal(0, sigma_var ** 2)) * math.exp(np.random.normal(0, sigma_var ** 2))
                        mutated_ind['dim'][key] = parent['dim'][key] + np.random.normal(0, mutated_ind['sigma'][key])

                    mutated_ind = self.validate(mutated_ind)
                    mutated_ind['eval'] = self.evaluate(self.evaluation_expression,
                                                        **mutated_ind['dim'])  # faz avaliação
                    offspring.append(mutated_ind)

            offspring_and_parents = population + offspring
            best = sorted(offspring_and_parents, key=lambda i: i['eval'])
            population = best[0:num_parents]

        return population

    def one_plus_one_ES_test(self, sigma=1.0, c=0.817, n=10, iter=100, seed=0):
        self.evaluation_counter = 0
        np.random.seed(seed)
        p = []
        A = []
        t = 0
        ind = {'dim': [0] * self.num_dimensions}
        k_checkpoint = 0
        self.reset_execution_history()

        for d in range(self.num_dimensions):
            ind['dim'][d] = np.random.uniform(self.get_bound_value_test(d, "lower"), self.get_bound_value_test(d, "upper"), 1)[0]

        ind = self.validate_test(ind)
        ind['eval'] = self.evaluate_test(ind)
        p.append(ind)

        for i in range(iter):
            t += 1
            mutated_ind = {'dim': [0] * self.num_dimensions}

            # mutate individual
            for d in range(self.num_dimensions):
                mutated_ind['dim'][d] = ind['dim'][d] + sigma * np.random.normal()

            mutated_ind = self.validate_test(mutated_ind)
            # calculate new phi based on the mutation
            mutated_ind['eval'] = self.evaluate_test(mutated_ind)

            if mutated_ind['eval'] < ind['eval']:
                ind = mutated_ind
                A.append(1)
                p.append(mutated_ind)
            else:
                A.append(0)

            if t % n == 0:
                # get successes and failures from at most 10n entries in A
                window = A[-10 * n:]
                success = sum(window)
                failure = abs(success - len(window))
                ps = success / (success + failure)
                if ps < 1 / 5:
                    sigma = sigma * c
                elif ps > 1 / 5:
                    sigma = sigma / c

            if self.get_evaluations_count() >= self.calculate_checkpoint(k_checkpoint, iter):
                error = self.abs_error(p[len(p) - 1]['eval'])
                self.add_checkpoint(self.get_evaluations_count(), error)
                k_checkpoint += 1

            if len(p)>=2 and  self.stop(p[len(p) - 1]['eval']):
                break

        self.add_checkpoint(self.get_evaluations_count(), self.abs_error(p[len(p) - 1]['eval']))
        return p

    def populational_isotropic_ES_test(self, dimension_gen_interval=(0, 0), sigma_var=0.5, iter=100, seed=0,
                                       num_parents=0, num_offspring=0, sigma_interval=(0.1, 10)):
        self.evaluation_counter = 0
        np.random.seed(seed)
        population = []
        t = 0
        k_checkpoint = 0
        self.reset_execution_history()

        for i in range(num_parents):
            individual = {'dim': [0] * self.num_dimensions, 'sigma': None}
            for d in range(self.num_dimensions):
                # individual['dim'][d] = np.random.uniform(dimension_gen_interval[0], dimension_gen_interval[1], 1)[0]
                individual['dim'][d] = np.random.uniform(self.get_bound_value_test(d, "lower"), self.get_bound_value_test(d, "upper"), 1)[0]

            # individual['sigma'] = np.random.uniform(0, 1, 1)[0]
            individual['sigma'] = np.random.uniform(sigma_interval[0], sigma_interval[1], 1)[0]

            individual = self.validate_test(individual)
            evaluation = self.evaluate_test(individual)
            individual['eval'] = evaluation
            population.append(individual)

        # todo: mudar criterio para chamadas da função objetivo talvez

        for i in range(iter):
            t += 1
            offspring = []

            for parent in population:
                for j in range(math.ceil(num_offspring / num_parents)):
                    mutated_ind = {'dim': [0] * self.num_dimensions, 'sigma': None}
                    mutated_ind['sigma'] = parent['sigma'] * math.exp(
                        np.random.normal(0, sigma_var ** 2))  # muta o sigma
                    for d in range(self.num_dimensions):  # muta cada uma das dimensões e seus sigmas
                        mutated_ind['dim'][d] = parent['dim'][d] + np.random.normal(0, mutated_ind['sigma'])
                    mutated_ind = self.validate_test(mutated_ind)
                    mutated_ind['eval'] = self.evaluate_test(mutated_ind)  # faz avaliação
                    offspring.append(mutated_ind)

            offspring_and_parents = population + offspring
            best = sorted(offspring_and_parents, key=lambda i: i['eval'])
            population = best[0:num_parents]

            if self.get_evaluations_count() >= self.calculate_checkpoint(k_checkpoint, iter*num_offspring):
                error = self.abs_error(population[0]['eval'])
                self.add_checkpoint(self.get_evaluations_count(), error)
                k_checkpoint += 1

            if self.stop(population[0]['eval']):
                break

        self.add_checkpoint(self.get_evaluations_count(), self.abs_error(population[0]['eval']))
        return population

    def populational_non_isotropic_ES_test(self, dimension_gen_interval=(0, 0), sigma_var=0.5, iter=100, seed=0,
                                           num_parents=0, num_offspring=0, sigma_interval=(0.1, 10)):
        self.evaluation_counter = 0
        np.random.seed(seed)
        population = []
        t = 0
        k_checkpoint = 0
        self.reset_execution_history()

        for i in range(num_parents):
            individual = {'dim': [0] * self.num_dimensions, 'sigma': [0] * self.num_dimensions}
            for d in range(self.num_dimensions):
                # individual['dim'][d] = np.random.uniform(dimension_gen_interval[0], dimension_gen_interval[1], 1)[0]
                individual['dim'][d] = np.random.uniform(self.get_bound_value_test(d, "lower"), self.get_bound_value_test(d, "upper"), 1)[0]
                # individual['sigma'][d] = np.random.uniform(0, 1, 1)[0]
                individual['sigma'][d] = np.random.uniform(sigma_interval[0], sigma_interval[1], 1)[0]

            individual = self.validate_test(individual)
            evaluation = self.evaluate_test(individual)
            individual['eval'] = evaluation
            population.append(individual)

        # todo: mudar criterio para chamadas da função objetivo talvez

        for i in range(iter):
            t += 1
            offspring = []

            for parent in population:
                for j in range(math.ceil(num_offspring / num_parents)):
                    mutated_ind = {'dim': [0] * self.num_dimensions, 'sigma': [0] * self.num_dimensions}
                    for i in range(len(parent['dim'])):  # muta cada uma das dimensões e seus sigmas
                        mutated_ind['sigma'][i] = parent['sigma'][i] * math.exp(
                            np.random.normal(0, sigma_var ** 2)) * math.exp(np.random.normal(0, sigma_var ** 2))
                        mutated_ind['dim'][i] = parent['dim'][i] + np.random.normal(0, mutated_ind['sigma'][i])

                    mutated_ind = self.validate_test(mutated_ind)
                    mutated_ind['eval'] = self.evaluate_test(mutated_ind)  # faz avaliação
                    offspring.append(mutated_ind)

            offspring_and_parents = population + offspring
            best = sorted(offspring_and_parents, key=lambda i: i['eval'])
            population = best[0:num_parents]

            if self.get_evaluations_count() >= self.calculate_checkpoint(k_checkpoint, iter*num_offspring):
                error = self.abs_error(population[0]['eval'])
                self.add_checkpoint(self.get_evaluations_count(), error)
                k_checkpoint += 1

            if self.stop(population[0]['eval']):
                break

        self.add_checkpoint(self.get_evaluations_count(), self.abs_error(population[0]['eval']))
        return population
