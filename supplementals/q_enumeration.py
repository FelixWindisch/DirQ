import numpy as np
import random
from itertools import combinations, chain
import seaborn as sns
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, Normalize

def generate_random_unidirectional_simplex_pair(n, m, q):
    vertex_set_1 = [i for i in range(n)]
    vertex_set_2 = [i for i in range(m)]
    return (random.sample(vertex_set_1, q+1), random.sample(vertex_set_2, q+1))


def enumerate_unidirectional_simplex_pairs(n, m, q, shared_face_dimension):
    vertex_set_1 = [i for i in range(n)]
    vertex_set_2 = [i for i in range(m)]
    result = []

    for c1 in combinations(vertex_set_1, shared_face_dimension):
        for c2 in combinations(vertex_set_2, shared_face_dimension):
            result.append((list(c1), list(c2)))
    return result
def is_q_old_near(s1, s2, n, m, q, i, j):
    if len(s1) > q+2:
        return True
    #if i > len(s1):
    #    i = n-1
    #if j > len(s2):
    #    j = m-1
    if len(s1) == q+2:
        if i in s1 and j in s2:
            return False
    else:
        if i in s1 or j in s2:
            return False
    return True


def is_q_new_near(s1, s2, n, m, q, i, j):
    if len(s1) > q+2:
        return True
    if len(s1) == q+2:
        if i == 0:
            s1.append(-1)
        if j == 0:
            s2.append(-1)
        if i >= len(s1):
            i = len(s1)
            s1.append(n)
        if j >= len(s2):
            j = len(s2)
            s2.append(m)
        if s1[i] > (s1[i - 1] + 1) or s2[j] > (s2[j - 1] + 1):
            return True
    if i == 0:
        s1.append(-1)
    if j == 0:
        s2.append(-1)
    if i >= len(s1):
        i = len(s1)
        s1.append(n)
    if j >= len(s2):
        j = len(s2)
        s2.append(m)


    if s1[i] > (s1[i-1]+1) and s2[j] > (s2[j-1]+1):
        return True
    return False
def plot_ij_both(n, m, q, shared_face_dimension):
    results_new = np.zeros((n, m))
    results_old = np.zeros((n, m))
    results_both = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            enumeration = enumerate_unidirectional_simplex_pairs(n, m, q, shared_face_dimension)
            for (s1, s2) in enumeration:
                #s1, s2 = generate_random_unidirectional_simplex_pair(n, m, q)
                s1.sort()
                s2.sort()
                q_old_near = is_q_old_near(s1, s2, n, m, q, i, j)
                q_new_near = is_q_new_near(s1, s2, n, m, q, i, j)
                if q_old_near:
                    results_old[i][j] += 1
                if q_new_near:
                    results_new[i][j] += 1
                if q_new_near and q_old_near:
                    results_both[i][j] += 1
    results_old /= len(enumeration)
    results_new /= len(enumeration)
    results_both /= len(enumeration)

    results_both_norm = results_both - results_both.min()
    results_both_norm = results_both_norm / results_both_norm.max()
    sns.color_palette("crest", as_cmap=True)
    ax = sns.heatmap(results_both_norm, linewidth=0.5, cmap="crest")
    ax.set_xlabel('i', fontsize=20)
    ax.set_ylabel('j', fontsize=20)
    ax.invert_yaxis()

    plt.show()

    print(results_old)
    print(results_new)
    print(results_both)


def plot_by_size(n_min, n_max, q, shared_face_dimension):
    results_new = np.zeros((n_max-n_min, n_max-n_min))
    results_old = np.zeros((n_max-n_min, n_max-n_min))
    for n in range(n_min, n_max):
        for m in range(n_min, n_max):
            enumeration = enumerate_unidirectional_simplex_pairs(n, m, q, shared_face_dimension)
            for (s1, s2) in enumeration:
                #s1, s2 = generate_random_unidirectional_simplex_pair(n, m, q)
                s1.sort()
                s2.sort()
                q_old_near = is_q_old_near(s1, s2, n, m, q, 0, 0)
                q_new_near = is_q_new_near(s1, s2, n, m, q, 0, 0)
                if q_old_near:
                    results_old[n-n_min][m-n_min] += 1
                if q_new_near:
                    results_new[n-n_min][m-n_min] += 1
            results_old[n - n_min][m - n_min] /= len(enumeration)
    #results_old /= (n_max-n_min) * (n_max-n_min)
    #results_new /= (n_max-n_min) * (n_max-n_min)

    results_both_norm = results_old - results_old.min()
    results_both_norm = results_both_norm / results_both_norm.max()
    sns.color_palette("crest", as_cmap=True)
    x_axis_labels = list(range(n_min, n_max))
    y_axis_labels = list(range(n_min, n_max))
    ax = sns.heatmap(results_old, linewidth=0.5, cmap="crest", xticklabels=x_axis_labels, yticklabels=y_axis_labels)#, norm=LogNorm())
    ax.set_xlabel('dim(τ)', fontsize=14)
    ax.set_ylabel('dim(σ)', fontsize=14)
    ax.invert_yaxis()

    plt.show()

q = 4
plot_by_size(6, 14, q, q+1)
plot_by_size(6, 14, q, q+2)
q = 3
plot_ij_both(9,9, q, q+1)
plot_ij_both(9,9, q, q+2)
q = 5
plot_ij_both(9,9, q, q+1)
plot_ij_both(9,9, q, q+2)

