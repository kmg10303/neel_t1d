Loading data
Fetching IL-6 exposure data for 5 instruments...
Found 5 genetic instruments.

----- Instrument Strength (F-statistics) -----
     rsid  beta_exp   se_exp    F_stat  weak
rs7553796  0.000893 0.001856  0.231498  True
rs4537545 -0.001134 0.001855  0.373713  True
rs2228145 -0.001057 0.001868  0.320181  True
rs4329505  0.001859 0.002496  0.554715  True
rs1800795 -0.015203 0.001863 66.593660 False
WARNING: 4 instrument(s) have F < 10 (weak instrument bias risk).

----- MR Results (IVW via WLS) -----
/Users/kavimongia-gasper/neel_research_index/venv/lib/python3.9/site-packages/statsmodels/stats/stattools.py:74: ValueWarning: omni_normtest is not valid with less than 8 observations; 5 samples were given.
  warn("omni_normtest is not valid with less than 8 observations; %i "
                                 WLS Regression Results                                
=======================================================================================
Dep. Variable:               beta_out   R-squared (uncentered):                   0.247
Model:                            WLS   Adj. R-squared (uncentered):              0.058
Method:                 Least Squares   F-statistic:                              1.311
Date:                Sun, 31 May 2026   Prob (F-statistic):                       0.316
Time:                        16:11:30   Log-Likelihood:                          9.6038
No. Observations:                   5   AIC:                                     -17.21
Df Residuals:                       4   BIC:                                     -17.60
Df Model:                           1                                                  
Covariance Type:            nonrobust                                                  
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
beta_exp       2.7490      2.401      1.145      0.316      -3.918       9.416
==============================================================================
Omnibus:                          nan   Durbin-Watson:                   0.167
Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.548
Skew:                           0.611   Prob(JB):                        0.760
Kurtosis:                       1.934   Cond. No.                         1.00
==============================================================================

Notes:
[1] R² is computed without centering (uncentered) since the model does not contain a constant.
[2] Standard Errors assume that the covariance matrix of the errors is correctly specified.
Done