# Complexity.sage

## Input

```python
TOTALTIME(n,k,m,q,sigma,beta_range,itr_range,P,qc,mlwe,qset='s')
```

* (n,k,m,q,sigma) : Module-LWE parameter
    * lattice size : $n*k \times n*m$, $d=n*m$

* $beta$ : BKZ blocksize
* itr : repetition number of BKZ

* P : total success probability 

* qc：
    * q : Quantum
    * c : Classical

* mlwe : 
    * n : without algebraic structures
    * m : with algebraic structures

* qset：
    * p : quantum multi-target BDD (Proposed method)
    * s : naive method
    


## Output
itr,beta: value that minimize TOTAL 

* Classical : [itr, beta, TOTAL, BKZ, ENUM]
* Quantum : [itr, beta, TOTAL, BKZ, DETECT, FIND]


