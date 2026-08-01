cause they are not in the dataset. Configure the logging level to `logging.WARNING` or higher for additional details.
  warnings.warn(
STEP 2 - Identification (what is estimable given the DAG?)
Estimand type: EstimandType.NONPARAMETRIC_ATE

### Estimand : 1
Estimand name: backdoor
No such variable(s) found!

### Estimand : 2
Estimand name: iv
Estimand expression:
 ⎡                                                                                   -1⎤
 ⎢                d                       ⎛                d                        ⎞  ⎥
E⎢──────────────────────────────────(T1D)⋅⎜──────────────────────────────────([IL₆])⎟  ⎥
 ⎣d[rs₂₀₆₉₈₅₂  rs₁₅₂₄₁₀₇  rs₁₅₅₄₆₀₆]      ⎝d[rs₂₀₆₉₈₅₂  rs₁₅₂₄₁₀₇  rs₁₅₅₄₆₀₆]       ⎠  ⎦
Estimand assumption 1, As-if-random: If U→→T1D then ¬(U →→{rs2069852,rs1524107,rs1554606})
Estimand assumption 2, Exclusion: If we remove {rs2069852,rs1524107,rs1554606}→{IL6}, then ¬({rs2069852,rs1524107,rs1554606}→T1D)

### Estimand : 3
Estimand name: frontdoor
No such variable(s) found!


  Naive regression (adjusts BMI/AGE/SMOKING only) : 4.1582
  DoWhy IV estimate                              : 3.2544
  IVW/WLS from mr_t1d.py (summary-level MR)       : 3.4829 (95% CI 2.107 to 4.859)

  Naive bias vs truth : +0.6753
  DoWhy IV bias vs truth : -0.2285
  -> IV estimate lands inside the IVW confidence interval. The DoWhy
     IV estimator and our IVW regression target the same estimand.