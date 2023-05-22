from ES import ESAlgorithm
import time
from py_expression_eval import Parser
from GRN_aux_5 import evaluation_5
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
    alg = ESAlgorithm()
    alg.set_evaluation_function(evaluation_5)
    alg.set_num_dimensions(IND_SIZE)
    best_inds = [None, None, None, None, None]

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

    res_1 = alg.one_plus_one_ES_test(sigma=0.5, c=0.817, n=10, iter=20000, seed=seed+1)
    if best_inds[0] is None or res_1[-2:-1][0]['eval'] < best_inds[0][-2:-1][0]['eval']:
        best_inds[0] = res_1

    res_2 = alg.populational_isotropic_ES_test(sigma_var=0.5, dimension_gen_interval=(-10, 10), iter=1000, seed=seed+2,num_parents=10, num_offspring=20)
    if best_inds[1] is None or res_2[0]['eval'] < best_inds[1][0]['eval']:
        best_inds[1] = res_2

    res_3 = alg.populational_non_isotropic_ES_test(sigma_var=0.5, dimension_gen_interval=(-10, 10), iter=1000, seed=seed+3, num_parents=10, num_offspring=20)
    if best_inds[2] is None or res_3[0]['eval'] < best_inds[2][0]['eval']:
        best_inds[2] = res_3

    res_4 = alg.populational_isotropic_ES_test(sigma_var=0.5, dimension_gen_interval=(-10, 10), iter=2000, seed=seed+4,num_parents=5, num_offspring=10)
    if best_inds[3] is None or res_4[0]['eval'] < best_inds[3][0]['eval']:
        best_inds[3] = res_4

    res_5 = alg.populational_non_isotropic_ES_test(sigma_var=0.5, dimension_gen_interval=(-10, 10), iter=2000, seed=seed+5, num_parents=5, num_offspring=10)
    if best_inds[4] is None or res_5[0]['eval'] < best_inds[4][0]['eval']:
        best_inds[4] = res_5

    res_str = f"{ res_1[-2:-1][0]['eval']},{res_2[0]['eval']},{res_3[0]['eval']},{res_4[0]['eval']},{res_5[0]['eval']}"

    return res_str, best_inds


def run_experiment():
    num_runs = 30
    result_file = open("../results/exp1/GRN5-20000-python_impl.csv", "a")
    time_file = open("../results/exp1/GRN5-20000-python_impl-time.csv", "a")
    best_inds_file = open("../results/exp1/GRN5-20000-python_impl-best_inds.txt", "a")
    result_file.write("1+1,10+20-i,10+20-ni,5+10-i,5+10-ni\n")

    for i in range(num_runs):
        print(f"Run {i}")

        start = time.time()
        exp_result, best_inds = round(i)
        csv_str_results = f"{exp_result}\n"
        end = time.time()

        print(f"Time elapsed: {end - start}")
        csv_str_time = f"{end - start}\n"

        result_file.write(csv_str_results)
        result_file.flush()
        time_file.write(csv_str_time)
        time_file.flush()
        best_inds_file.write(str(best_inds) + "\n\n")
        best_inds_file.flush()

    result_file.close()
    time_file.close()
    best_inds_file.close()


if __name__ == '__main__':
    run_experiment()
    #result_file = open("../results/exp1/testando.txt", "a")
    #result_file.write("testando")
    #result_file.close()
