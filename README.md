## Replicating Arcidiacono Miller (2011, Table 1). 

The code is structured to be legible rather than efficient. Notebooks are meant to be self-contained. 

The code for FIML estimators is provided, but I do not run the Monte Carlos given the long expected computation time. 

#### Notebooks
- ```ArcidiaconoMiller_Case1```: _To obtain the estimates reported in column 4 of Table I (when s is ignored), we estimate the CCP’s using_ $W_{1t} ≡ (1, x_{1t}, x^2_{1t}, x_2, x^2_2, x_{1t} x_2) _as regressors in a logit._ (Supp. Mat, p.12)
- ```ArcidiaconoMiller_Case2```: _For the parameters reported in column 3, $W_{1t}$ is fully interacted with $W_{2t} ≡ (1, s, t, st, t^2, st^2)$, which is 36 parameters to estimate in the logit generating the CCP’s. Since s is observed, this flexible logit is estimated once._ (Supp. Mat, p.13)

#### Content 
```7743_data and programs_0``` is the replication material downloaded from the Econometrica webpage on April, 1st 2025.

In the replication material, I believe
- ```nohetero``` refers to column (3) of Table 1 on page 1854 [Case 2 in Section B.1.3 of the Supp Mat]
- ```noheterowrong``` refers to column (4) of Table 1 on page 1854
- ```hetero10ns``` refers to column (8) of Table 1 on page 1854 (note absence of FIML in folder)
- ```nohetero10ns``` refers to column (7) of Table 1 on page 1854 (note absence of FIML in folder)

#### Details 
Notebooks are tested on Python 3.12.9. Times based on running the Jupyter notebook in VS Code on Windows 10. 
Processor: AMD Ryzen 3 3100 4-Core Processor 3.59 GHz
Installed RAM	16.0 GB