from ES import populational_isotropic_ES, populational_non_isotropic_ES, ESAlgorithm
import matplotlib.pyplot as plt
import numpy as np
import time
from py_expression_eval import Parser
from GRN.GRN_aux_10 import evaluation
from GRN.GRN_aux_5 import evaluation_5
parser = Parser()

IND_SIZE    = 15  # Tamanho do indivíduo (quantidade de coeficientes)
MIN_K       = 0.1  # Menor valor que K pode assumir
MAX_K       = 1  # Maior valor que K pode assumir
MIN_N       = 1  # Menor valor que N pode assumir
MAX_N       = 25  # Maior valor que N pode assumir
MIN_TAU     = 0.1  # Menor valor que TAU pode assumir
MAX_TAU     = 5  # Maior valor que TAU pode assumir
MIN_STRATEGY = 0.1  # Menor valor que a estratégia pode assumir
MAX_STRATEGY = 10  # Maior valor que a estratégia pode assumir
TAU_SIZE    = 5
N_SIZE      = 5
K_SIZE      = 5

def round(seed=0):
    experiment = {}
    error_history = {'1p1': None, 'pi': None, 'pni': None}

    alg = ESAlgorithm()
    alg.set_evaluation_function(evaluation_5)
    alg.set_num_dimensions(IND_SIZE)
    #alg.set_global_variable_bounds(100, True, -100, True)

    cont = 0
    for j in range(0, TAU_SIZE, 1):
        alg.set_variable_bounds_test(j, MAX_TAU, True, MIN_TAU, True)
        cont = j

    for j in range(cont + 1, cont + 1 + K_SIZE, 1):
        alg.set_variable_bounds_test(j, MAX_K, True, MIN_K, True)
        cont = j

    for j in range(cont + 1, cont + 1 + N_SIZE, 1):
        alg.set_variable_bounds_test(j, MAX_N, True, MIN_N, True)
        cont = j

    #alg.set_error_stop_criterion(0.0001)

    res_1 = alg.one_plus_one_ES_test(sigma=0.5, c=0.817, n=10, iter=2000, seed=seed)
    error_history['1p1'] = alg.get_execution_history()

    res_2 = alg.populational_isotropic_ES_test(sigma_var=0.5, dimension_gen_interval=(-10, 10), iter=100, seed=seed,num_parents=10, num_offspring=20)
    error_history['pi'] = alg.get_execution_history()

    res_3 = alg.populational_non_isotropic_ES_test(sigma_var=0.5, dimension_gen_interval=(-10, 10), iter=100, seed=seed, num_parents=10, num_offspring=20)

    error_history['pni'] = alg.get_execution_history()

    experiment['1p1'] = res_1[-2:-1][0]['eval']
    experiment['pi'] = res_2[0]['eval']
    experiment['pni'] = res_3[0]['eval']

    res_str = f"{ res_1[-2:-1][0]['eval']},{res_2[0]['eval']},{res_3[0]['eval']}"

    return res_str


def run_experiment():
    num_exp = 3
    result_file = open("../results/GRN5-2000-python_impl.csv", "a")
    time_file = open("../results/GRN5-2000-python_impl-time.csv", "a")
    result_file.write("1p1,pi,pni\n")

    for i in range(num_exp):
        print(f"Run {i}")

        start = time.time()
        csv_str_results = f"{round(i)}\n"
        end = time.time()

        print(f"Time elapsed: {end - start}")
        csv_str_time = f"{end - start}\n"

        result_file.write(csv_str_results)
        time_file.write(csv_str_time)

    result_file.close()
    time_file.close()

if __name__ == '__main__':
    run_experiment()