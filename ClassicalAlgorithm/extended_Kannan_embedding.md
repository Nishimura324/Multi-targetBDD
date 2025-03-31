# Extended Kannan Embedding

## Input

```python
sol_Kannan(a, t, q, sigma, eta, r, beta, loop_max)
```

* $(a, t)$ : Module-LWE instance

* $(q,\sigma)$ : Module-LWE parameter

* $\eta$ : lower right diagonal component of embedding matrix

* $r$ : the number of embedding vector

* $\beta$ : BKZ blocksize

* $\mathtt{loop} \ \mathtt{max}$ : BKZ tour steps

## Output

```python
return e
```

* $e$ : Module-LWE error polynomial（Output 0 if not find

* $\displaystyle e_i(x) = \sum_{\ell=1}^{n} e\lbrack i\rbrack \lbrack \ell\rbrack\times x^\ell$

## Note

* In BKZ, We did not dare to introduce AUTO_ABORT. BKZ option reference is  <https://www.maths.ox.ac.uk/system/files/attachments/lab-02.pdf> .

* When parameter $nm$ is large (over 200?), LLL precision=double is insufficient, so

    ```python
    B = BKZ.reduction(B, par, float_type = "ld")
    ```

    But Run time is longer.
