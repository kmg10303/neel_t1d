----- Instrument Strength (F-statistics) -----
     rsid  beta_exp   se_exp    F_stat  weak
rs1800797 -0.015203 0.001863 66.593660 False
rs1800795 -0.015203 0.001863 66.593660 False
rs2069837 -0.001624 0.003649  0.198073  True
rs1554606 -0.015253 0.001862 67.104429 False
rs1524107 -0.025149 0.004419 32.388694 False
WARNING: 1 instrument(s) have F < 10 (weak instrument bias risk).

----- MR Results (IVW via WLS) -----
/Users/kavimongia-gasper/neel_research_index/venv/lib/python3.9/site-packages/statsmodels/stats/stattools.py:74: ValueWarning: omni_normtest is not valid with less than 8 observations; 5 samples were given.
  warn("omni_normtest is not valid with less than 8 observations; %i "
                                 WLS Regression Results                                
=======================================================================================
Dep. Variable:               beta_out   R-squared (uncentered):                   0.885
Model:                            WLS   Adj. R-squared (uncentered):              0.856
Method:                 Least Squares   F-statistic:                              30.71
Date:                Sun, 31 May 2026   Prob (F-statistic):                     0.00518
Time:                        16:36:13   Log-Likelihood:                          12.810
No. Observations:                   5   AIC:                                     -23.62
Df Residuals:                       4   BIC:                                     -24.01
Df Model:                           1                                                  
Covariance Type:            nonrobust                                                  
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
beta_exp       2.9792      0.538      5.542      0.005       1.487       4.472
==============================================================================
Omnibus:                          nan   Durbin-Watson:                   2.226
Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.826
Skew:                          -0.415   Prob(JB):                        0.662
Kurtosis:                       1.190   Cond. No.                         1.00
==============================================================================

Notes:
[1] R² is computed without centering (uncentered) since the model does not contain a constant.
[2] Standard Errors assume that the covariance matrix of the errors is correctly specified.
Done