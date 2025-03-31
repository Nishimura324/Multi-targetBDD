from fpylll import *



def sol_Kannan(a, t, q, sigma, eta, r, beta, loop_max):
    m, k, n = len(a), len(a[0]), len(a[0][0])
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
    B = IntegerMatrix(d + r, d + r)
    for x in range(d):
        for y in range(d):
            B[x, y] = A[n * k + x, y]

    for h in range(r):
        for i in range(m):
            for l in range(n):
                if h <= l:
                    B[d + h, n * i + l] = t[i][l - h]
                else:
                    B[d + h, n * i + l] = -t[i][l - h + n]

        B[d + h, d + h] = eta

    flags = BKZ.MAX_LOOPS | BKZ.GH_BND
    par = BKZ.Param(beta, strategies=BKZ.DEFAULT_STRATEGY,
                    max_loops=loop_max, flags=flags)
    B = BKZ.reduction(B, par)

    e = [[0 for l in range(n)] for i in range(m)]

    for h in range(r):
        flag = True
        for i in range(r):
            if i != h and B[0][d + i] != 0:
                flag = False

        if abs(B[0][d + h]) != eta:
            flag = False

        sgn = B[0][d + h] // eta

        if flag:
            for i in range(m):
                for l in range(n):
                    if l + h < n:
                        e[i][l] = B[0][n * i + l + h] * sgn
                    else:
                        e[i][l] = -B[0][n * i + l + h - n] * sgn

            return e

    return e
