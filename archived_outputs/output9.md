Loading data
Fetching IL-6 exposure data for 5 instruments...
Found 5 genetic instruments.

----- Instrument Strength (F-statistics) -----
     rsid  beta_exp   se_exp    F_stat  weak
rs7802307 -0.015203 0.001863 66.593660 False
rs1554606 -0.015253 0.001862 67.104429 False
rs1800796 -0.025037 0.004394 32.467173 False
rs1524107 -0.025149 0.004419 32.388694 False
rs2069852 -0.026594 0.005321 24.979331 False

----- MR Results (IVW via WLS) -----
/Users/kavimongia-gasper/neel_research_index/venv/lib/python3.9/site-packages/statsmodels/stats/stattools.py:74: ValueWarning: omni_normtest is not valid with less than 8 observations; 5 samples were given.
  warn("omni_normtest is not valid with less than 8 observations; %i "
                                 WLS Regression Results                                
=======================================================================================
Dep. Variable:               beta_out   R-squared (uncentered):                   0.920
Model:                            WLS   Adj. R-squared (uncentered):              0.900
Method:                 Least Squares   F-statistic:                              45.90
Date:                Sun, 31 May 2026   Prob (F-statistic):                     0.00248
Time:                        16:54:30   Log-Likelihood:                          11.955
No. Observations:                   5   AIC:                                     -21.91
Df Residuals:                       4   BIC:                                     -22.30
Df Model:                           1                                                  
Covariance Type:            nonrobust                                                  
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
beta_exp       3.4414      0.508      6.775      0.002       2.031       4.852
==============================================================================
Omnibus:                          nan   Durbin-Watson:                   1.088
Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.693
Skew:                           0.226   Prob(JB):                        0.707
Kurtosis:                       1.233   Cond. No.                         1.00
==============================================================================

Notes:
[1] R² is computed without centering (uncentered) since the model does not contain a constant.
[2] Standard Errors assume that the covariance matrix of the errors is correctly specified.
Done