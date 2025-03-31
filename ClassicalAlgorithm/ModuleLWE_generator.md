# Module-LWE instance generator

## Input

```python
gen(n, k, m, q, sigma)
```

* $(n,k,m,q,\sigma)$ : Module-LWE parameter

## Output

```python
return a, s, e, t
```

* $(a,t)$ : Module-LWE instance

* $\displaystyle a_{i,j}(x) = \sum_{\ell=1}^{n} a\lbrack i\rbrack \lbrack j\rbrack\lbrack \ell\rbrack\times x^\ell$

* $\displaystyle s_j(x) = \sum_{\ell=1}^{n} s\lbrack j\rbrack \lbrack \ell\rbrack\times x^\ell$

* $\displaystyle e_i(x) = \sum_{\ell=1}^{n} e\lbrack i\rbrack \lbrack \ell\rbrack\times x^\ell$

* $\displaystyle t_i(x) = \sum_{\ell=1}^{n} t\lbrack i\rbrack \lbrack \ell\rbrack\times x^\ell$

## Method

* $(a, s)$ is generated randomly

* each $e$ coefficient is sampled from Discrete Gaussian Distribution$D_{\mathbb{Z},\sigma}$

* $\displaystyle t_i(x) = \sum_{j=1}^{k}a_{i,j}(x)s_j(x) + e_i(x)$
