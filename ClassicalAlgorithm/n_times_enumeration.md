# Proposed method

## Input

``` python
sol_enum(a, t, q, sigma, beta, loop_max, beta_prime, p0, T=300, t0=0.01, t1=0.0001, s=20000, r=-1)
```

* $(a, t)$ : Module-LWE instance

* $(q,\sigma)$ : Module-LWE parameter

* $\beta$ : BKZ blocksize

* $\mathtt{loop} \ \mathtt{max}$ : BKZ tour steps

* $\beta'$ : block size executing enumeration

* $p_0$ : target success probability

* $T$ : the number of error vector samples 

* $(t_0, t_1)$ : initial and final temperature in Simulated Annealing


* $s$ : Simulated Annealing steps

* $r$ : the number of target vector (basically $r = n$)

## Output

``` python
return e
```

* $e$ : Module-LWE error polynomial（Output 0 if not find）

* $\displaystyle e_i(x) = \sum_{\ell=1}^{n} e\lbrack i\rbrack \lbrack \ell\rbrack\times x^\ell$

## Preparation

* compile　[BDD_solver](../BDD_solver.cpp) using pybind11. Compile command is:

    ``` bash
    g++ -O3 -shared -std=c++17 -fPIC $(python3 -m pybind11 --includes) BDD_solver.cpp -o BDD_solver$(python3-config --extension-suffix)
    ```

    
