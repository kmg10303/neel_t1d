Loading data
Fetching IL-6 exposure data for 3 instruments...
Found 3 genetic instruments.

----- Instrument Strength (F-statistics) -----
     rsid  beta_exp   se_exp    F_stat  weak
rs1554606 -0.015253 0.001862 67.104429 False
rs1524107 -0.025149 0.004419 32.388694 False
rs2069852 -0.026594 0.005321 24.979331 False

----- MR Results (IVW via WLS) -----
/Users/kavimongia-gasper/neel_research_index/venv/lib/python3.9/site-packages/statsmodels/stats/stattools.py:74: ValueWarning: omni_normtest is not valid with less than 8 observations; 3 samples were given.
  warn("omni_normtest is not valid with less than 8 observations; %i "
                                 WLS Regression Results                                
=======================================================================================
Dep. Variable:               beta_out   R-squared (uncentered):                   0.925
Model:                            WLS   Adj. R-squared (uncentered):              0.887
Method:                 Least Squares   F-statistic:                              24.62
Date:                Sun, 31 May 2026   Prob (F-statistic):                      0.0383
Time:                        17:01:44   Log-Likelihood:                          7.1060
No. Observations:                   3   AIC:                                     -12.21
Df Residuals:                       2   BIC:                                     -13.11
Df Model:                           1                                                  
Covariance Type:            nonrobust                                                  
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
beta_exp       3.4829      0.702      4.962      0.038       0.463       6.503
==============================================================================
Omnibus:                          nan   Durbin-Watson:                   1.998
Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.376
Skew:                           0.435   Prob(JB):                        0.829
Kurtosis:                       1.500   Cond. No.                         1.00
==============================================================================

Notes:
[1] R² is computed without centering (uncentered) since the model does not contain a constant.
[2] Standard Errors assume that the covariance matrix of the errors is correctly specified.
Done