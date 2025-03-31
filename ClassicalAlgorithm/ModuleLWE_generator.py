import random
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

# a_{i,j}(x) = Σ[l=0,...,n-1] a[i][j][l]*x^l (1<=i<=m, 1<=j<=k)
# s_j(x) = Σ[l=0,...,n-1] s[j][l] * x ^ l (1<=j<=k)
# e_i(x) = Σ[l=0,...,n-1] e[i][l] * x ^ l (1<=i<=m)
# t_i(x) = Σ[l=0,...,n-1] t[i][l] * x ^ l (1<=i<=m)


def gen(n, k, m, q, sigma):
    a = [[[0 for l in range(n)] for j in range(k)] for i in range(m)]
    s = [[0 for l in range(n)] for j in range(k)]
    e = [[0 for l in range(n)] for i in range(m)]
    t = [[0 for l in range(n)] for i in range(m)]

    for i in range(m):
        for j in range(k):
            for l in range(n):
                a[i][j][l] = random.randint(-(q - 1) // 2, (q - 1) // 2)

    for j in range(k):
        for l in range(n):
            s[j][l] = random.randint(-(q - 1) // 2, (q - 1) // 2)

    for i in range(m):
        for l in range(n):
            e[i][l] = DiscreteGaussianDistributionIntegerSampler(sigma)()

    for i in range(m):
        for j in range(k):
            for l1 in range(n):
                for l2 in range(n):
                    if (l1 + l2 < n):
                        t[i][l1 + l2] += a[i][j][l1] * s[j][l2]
                    else:
                        t[i][l1 + l2 - n] -= a[i][j][l1] * s[j][l2]

        for l in range(n):
            t[i][l] += e[i][l]
            t[i][l] %= q
            if t[i][l] > (q - 1) // 2:
                t[i][l] -= q

    return a, s, e, t
