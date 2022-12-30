import matplotlib.pyplot as plt
import numpy as np


def plot_1(experiments):
    """
    Table of the worst, best, median, average, standard deviation of the experiment results for each algorithm.
    Also the boxplot of each algorithm's result set.
    :param experiments:
    :return:
    """
    res_1p1 = [e['1p1'] for e in experiments]
    res_pi = [e['pi'] for e in experiments]
    res_pni = [e['pni'] for e in experiments]
    tabela = [None, None, None]

    def table_row(array):
        row = [f"{np.average(array):.7f}", f"{np.min(array):.7f}", f"{np.max(array):.7f}", f"{np.median(array):.7f}", f"{np.std(array):.7f}"]
        return row

    tabela[0] = table_row(res_1p1)
    tabela[1] = table_row(res_pi)
    tabela[2] = table_row(res_pni)

    fig, ax = plt.subplots()
    # hide axes
    fig.patch.set_visible(False)

    ax.axis('off')
    ax.axis('tight')
    colLabels = ['Media', 'Melhor', 'Pior', 'Mediana', 'Desvio Padrão']
    rowLabels = ["(1+1)-ES", "Pop. Iso.", "Pop. N. Iso."]
    table = ax.table(tabela, colLabels=colLabels,rowLabels=rowLabels, loc='center')
    table.scale(3,2)
    table.set_fontsize(15)
    fig.tight_layout()

    fig_boxplot, ax_boxplot = plt.subplots()
    ax_boxplot.boxplot(x=[res_1p1, res_pi, res_pni], labels=["(1+1)-ES", "Pop. Iso.", "Pop. N. Iso."])

    plt.show()


def plot_2(error_history, title):
    """
    Plots the error history of the experiment for one of the algorithms
    :param experiment_error_history:
    :return:
    """
    #checkpoints = experiment_error_history['checkpoints']
    #values = experiment_error_history['values']

    fig, ax = plt.subplots()
    ax.plot(error_history)
    ax.set_title(title)
    plt.show()

