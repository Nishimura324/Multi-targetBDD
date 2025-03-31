from fpylll import *
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
import BDD_solver

# GSO matrix、GSO coefficients、GSO basis norm's square


def calc_GSO(B):
    d = len(B)
    GSO = [[0.0 for j in range(d)] for i in range(d)]
    mu = [[0.0 for j in range(d)] for i in range(d)]
    c = [0.0 for i in range(d)]

    for i in range(d):
        for j in range(d):
            GSO[i][j] += B[i][j]

    for i in range(d):
        for j in range(i):
            for k in range(d):
                mu[i][j] += GSO[i][k] * GSO[j][k]

            mu[i][j] /= c[j]
            for k in range(d):
                GSO[i][k] -= mu[i][j] * GSO[j][k]

        for j in range(d):
            c[i] += GSO[i][j] * GSO[i][j]

    return GSO, mu, c

# p0 : taget success probability
# T, t0, t1, s : Simulated annealing parameter in pruning parameter optimization
# enumerate for $r$ number of target vectors


def sol_enum(a, t, q, sigma, beta, loop_max, beta_prime, p0, T=300, t0=0.01, t1=0.0001, s=20000, r=-1):
    m, k, n = len(a), len(a[0]), len(a[0][0])
    if r == -1:
        r = n

    ############## reduction ##############

    d = n * m
    A = IntegerMatrix(n * k + d, d)
    for i in range(m):
        for j in range(k):
            for h in range(n):
                for l in range(n):
                    if (h <= l):
                        A[n * j + h, n * i + l] = a[i][j][l - h]
                    else:
                        A[n * j + h, n * i + l] = -a[i][j][l - h + n]

    for x in range(d):
        A[n * k + x, x] = q

    LLL.reduction(A)
    B = IntegerMatrix(d, d)
    for x in range(d):
        for y in range(d):
            B[x, y] = A[n * k + x, y]

    flags = BKZ.MAX_LOOPS | BKZ.GH_BND
    par = BKZ.Param(beta, strategies=BKZ.DEFAULT_STRATEGY,
                    max_loops=loop_max, flags=flags)
    B = BKZ.reduction(B, par)

    B_list = [[0 for y in range(d)] for x in range(d)]
    for x in range(d):
        for y in range(d):
            B_list[x][y] = B[x, y]

    GSO, mu, c = calc_GSO(B_list)

    ############## pruning ##############

    e_sample = [[0 for j in range(d)] for i in range(T)]
    for i in range(T):
        for j in range(d):
            e_sample[i][j] = DiscreteGaussianDistributionIntegerSampler(
                sigma)()

    R = BDD_solver.pruning_optimize(
        B_list, GSO, c, n, m, sigma, beta_prime, e_sample, p0, t0, t1, s, n)

    ############## enumeration ##############

    e = [[0 for l in range(n)] for i in range(m)]

    for h in range(r):
        rot_t = [0] * d
        for i in range(m):
            for l in range(n):
                if l + h < n:
                    rot_t[n * i + l + h] = t[i][l]
                else:
                    rot_t[n * i + l + h - n] = -t[i][l]

        rot_e = BDD_solver.sol(B_list, GSO, mu, c, rot_t, beta_prime, R)
        norm = 0
        for i in range(d):
            norm += rot_e[i] * rot_e[i]

        if norm != 0:
            for i in range(m):
                for l in range(n):
                    if l + h < n:
                        e[i][l] = rot_e[n * i + l + h]
                    else:
                        e[i][l] = -rot_e[n * i + l + h - n]

            return e

    return e
