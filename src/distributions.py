import random
import numpy as np
from scipy.stats import nbinom
from scipy.stats import powerlaw
from scipy.stats import pareto

def randomize_distr(distr, alpha: int = 2):
    sum_old = 0
    sum_new = 0
    for i in range(0, len(distr)):
        # weight = 1 / (i/2+1)
        # weight = 3/(i+1)
        weight = 0.3 + (10 / len(distr) * min(i, 10)) / 10
        random_element = random.uniform(-1 * weight, weight)
        distr[i] += abs(random_element * distr[i] * alpha)
        distr = list(distr)
        distr.sort(reverse=True)

    return distr

def linear_distr(data_points: int):
    return list(range(1, data_points + 1))

def pareto_distr(data_points: int):
    a = 0.15
    x_m = 1
    samples = np.linspace(start=0, stop=data_points + 1, num=data_points + 1)
    output = np.array(pareto.pdf(x=samples, b=a, loc=0, scale=x_m))

    values = randomize_distr(output.T[1:])
    scale = 1 / sum(values)

    values = [scale * val for val in values]

    return values
