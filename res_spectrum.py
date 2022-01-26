import scipy as sp
from scipy import special
import numpy as np
from typing import List
from sympy.utilities.iterables import multiset_permutations
import matplotlib.pyplot as plt
import pickle
import os.path as path
import os
import glob

def get_n_m(mu, m):
    return np.count_nonzero(mu == m)

def get_partitions(partitions_num, Lz, max_term):
    if Lz == 0:
        return []
    if partitions_num == 1:
        return [[i, Lz - i] for i in range(int(np.ceil(Lz / 2)), np.min([max_term + 1, Lz + 1]))]
    else:
        res = []
        for i in range(int(np.ceil(Lz / (partitions_num + 1))), np.min([max_term + 1, Lz + 1])):
            tails = get_partitions(partitions_num - 1, Lz - i, i - 1)
            res += [[i] + tail for tail in tails if i not in tail]
        return res


def sign(perm):
    N = len(perm)
    res = 1
    for i in range(N):
        id = np.where(np.array(perm) == i)[0][0]
        res *= (-1)**id
        perm.pop(id)
    return res


# TODO below is suitable for distinguishable particles, I think. So this may be a room for mistake.
def get_permutations(NA, NB, Lz, coefficients):
    partitions_num = NA + NB - 1
    partitions = [np.array(p) for p in get_partitions(partitions_num, Lz, Lz)]
    result = []
    for partition in partitions:
        for perm in multiset_permutations(list(range(NA+NB))):
            c = coefficients(partition[perm])
            if c != 0:
                result.append([partition[list(perm[:NA])], partition[list(perm[NA:])], sign(perm), c])
    return result


def get_F_B(m, N_phi):
    return special.betainc(m + 1, N_phi - m  + 1, np.cos(np.pi / 4)**2)


def get_R(NA: int, NB: int, Lz:int, coefficients, flux_quanta):
    permutations = get_permutations(NA, NB, Lz, coefficients)
    print('got permutations: ' + str(len(permutations)))
    mus = [list(map(int, item.split('-'))) for item in
           set([partition_list_to_str(permutations[i][0]) for i in range(len(permutations))])]
    nus = [list(map(int, item.split(' '))) for item in set([str(permutations[i][1]).replace('[', '').replace(']', '') for i in range(len(permutations))])]
    N_phi = flux_quanta(NA, NB)
    result = np.zeros((len(mus), len(nus)))
    # TODO change p_NA
    p_NA = np.math.factorial(NA + NB) / (np.math.factorial(NA) * np.math.factorial(NB)) / 2**(NA+NB)
    for permutation in permutations:
        mu = list(permutation[0])
        mu_i = mus.index(mu)
        nu = list(permutation[1])
        nu_i = nus.index(nu)
        if mu == [1, 0] and nu == [4, 2, 3]:
            b = 1
        factor = np.sqrt(
            np.prod([special.comb(get_n_m(mu, m) + get_n_m(nu, m), get_n_m(mu, m)) for m in range(N_phi)])) * \
                np.sqrt(np.prod([get_F_B(nu[i], N_phi) for i in range(NB)]))
        sign = permutation[2]
        c = permutation[3]
        result[mu_i, nu_i] += 1 / np.sqrt(p_NA) * c * sign * factor
    print('got R')
    return result, mus


def get_RA(NA: int, NB: int, Lz:int, coefficients, flux_quanta):
    R, mus = get_R(NA, NB, Lz, coefficients, flux_quanta)
    N_phi = flux_quanta(NA, NB)
    for mu_i in range(len(mus)):
        mu = mus[mu_i]
        factor = np.prod([np.sqrt(1 - get_F_B(mu[i], N_phi)) for i in range(NA)])
        R[mu_i, :] *= factor
    return R, mus


def get_rdm(NA: int, NB: int, Lz:int, coefficients, flux_quanta):
    RA, mus = get_RA(NA, NB, Lz, coefficients, flux_quanta)
    return np.matmul(RA, np.conj(np.transpose(RA))), mus


def iqh_1_coefficients(partition):
    if set(partition) == set(list(range(len(partition)))):
        return 1
    else:
        return 0

def iqh1_flux_quanta(NA, NB):
    return NA + NB - 1


def partition_list_to_str(partition):
    res = ''
    for item in partition:
        res += '-' + str(item)
    return res[1:]


def partition_str_to_list(partition):
    return [int(i) for i in partition.split('-')]

def differentiate(pairs, i):
    for pair_ind in range(len(pairs)):
        pair = pairs[pair_ind]
        if i in pair[:2]:
            pairs[pair_ind][3] *= pair[2] * (-1)**int(pair[1] == i)
            pairs[pair_ind][2] -= 1



def laughlin_3_coefficients(partition, nu=3):
    N = len(partition)



def laughlin_flux_quanta(NA, NB):
    N = NA + NB
    if NA == 4 and NB == 4:
        return 3

def initial(i, j, N, nu):
    return [[partition_list_to_str([0] * (i) + [3] + [0] * (N - 1 - i)), 1], \
            [partition_list_to_str([0] * (i) + [2] + [0] * (j - 1 - i) + [1] + [0] * (N-1-j)), -3], \
            [partition_list_to_str([0] * (i) + [1] + [0] * (j - 1 - i) + [2] + [0] * (N-1-j)), 3], \
            [partition_list_to_str([0] * (j) + [3] + [0] * (N - 1 - j)), -1]]

def multiply(partition1, partition2):
    p1 = partition_str_to_list(partition1[0])
    p2 = partition_str_to_list(partition2[0])
    return [partition_list_to_str([p1[i] + p2[i] for i in range(len(p1))]), partition1[1] * partition2[1]]

def laughlin_coefficients_precalc(N, nu=3):
    max_range = list(range(min([n for n in range(N+1) if (n+1)*n/2 >= N * 3]), nu*N+1))


def nu_1_es(NA, NB):
    N = NA + NB
    rho, mus = get_rdm(NA, NB, sum([i for i in range(N)]), iqh_1_coefficients, iqh1_flux_quanta)
    LzAs = np.array([sum(mu) for mu in mus])
    for l in range(1, 2 * N):
        inds = LzAs == l
        rho_l = rho[:, inds][inds]
        vals = np.linalg.eigvalsh(rho_l)
        for val in vals:
            plt.scatter(l - NA * (NA + NB - 1) / 2, -np.log(val), color='blue')
    plt.show()
nu_1_es(5, 6)

def laughlin_3_es(NA, NB):
    N = NA + NB
    rho, mus = get_rdm(NA, NB, sum([i for i in range(N)]), laughlin_3_coefficients, laughlin_flux_quanta)
    LzAs = np.array([sum(mu) for mu in mus])
    for l in range(1, 2 * N):
        inds = LzAs == l
        rho_l = rho[:, inds][inds]
        vals = np.linalg.eigvalsh(rho_l)
        for val in vals:
            plt.scatter(l - NA * (NA + NB - 1) / 2, -np.log(val), color='blue')
    plt.show()
# laughlin_3_es(4, 4)

