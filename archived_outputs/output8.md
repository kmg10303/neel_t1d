----- Instrument Strength (F-statistics) -----
     rsid  beta_exp   se_exp    F_stat  weak
rs1554606 -0.015253 0.001862 67.104429 False
rs1524107 -0.025149 0.004419 32.388694 False

----- MR Results (IVW via WLS) -----
/Users/kavimongia-gasper/neel_research_index/venv/lib/python3.9/site-packages/statsmodels/stats/stattools.py:74: ValueWarning: omni_normtest is not valid with less than 8 observations; 2 samples were given.
  warn("omni_normtest is not valid with less than 8 observations; %i "
                                 WLS Regression Results                                
=======================================================================================
Dep. Variable:               beta_out   R-squared (uncentered):                   0.910
Model:                            WLS   Adj. R-squared (uncentered):              0.819
Method:                 Least Squares   F-statistic:                              10.05
Date:                Sun, 31 May 2026   Prob (F-statistic):                       0.195
Time:                        16:50:28   Log-Likelihood:                          4.8578
No. Observations:                   2   AIC:                                     -7.716
Df Residuals:                       1   BIC:                                     -9.022
Df Model:                           1                                                  
Covariance Type:            nonrobust                                                  
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
beta_exp       3.3229      1.048      3.170      0.195      -9.994      16.640
==============================================================================
Omnibus:                          nan   Durbin-Watson:                   1.927
Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.333
Skew:                           0.000   Prob(JB):                        0.846
Kurtosis:                       1.000   Cond. No.                         1.00
==============================================================================

Notes:
[1] R² is computed without centering (uncentered) since the model does not contain a constant.
[2] Standard Errors assume that the covariance matrix of the errors is correctly specified.
Done