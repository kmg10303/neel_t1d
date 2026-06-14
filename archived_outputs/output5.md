Found 5 genetic instruments.

----- MR Results (IVW via WLS) -----
/Users/kavimongia-gasper/neel_research_index/venv/lib/python3.9/site-packages/statsmodels/stats/stattools.py:74: ValueWarning: omni_normtest is not valid with less than 8 observations; 5 samples were given.
  warn("omni_normtest is not valid with less than 8 observations; %i "
                                 WLS Regression Results                                
=======================================================================================
Dep. Variable:               beta_out   R-squared (uncentered):                   0.631
Model:                            WLS   Adj. R-squared (uncentered):              0.538
Method:                 Least Squares   F-statistic:                              6.826
Date:                Sun, 10 May 2026   Prob (F-statistic):                      0.0593
Time:                        16:49:33   Log-Likelihood:                         -7.2016
No. Observations:                   5   AIC:                                      16.40
Df Residuals:                       4   BIC:                                      16.01
Df Model:                           1                                                  
Covariance Type:            nonrobust                                                  
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
beta_exp      27.6771     10.593      2.613      0.059      -1.734      57.088
==============================================================================
Omnibus:                          nan   Durbin-Watson:                   1.821
Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.507
Skew:                          -0.067   Prob(JB):                        0.776
Kurtosis:                       1.446   Cond. No.                         1.00
==============================================================================

Notes:
[1] R² is computed without centering (uncentered) since the model does not contain a constant.
[2] Standard Errors assume that the covariance matrix of the errors is correctly specified.
/Users/kavimongia-gasper/neel_research_index/mr_t1d.py:136: FutureWarning: Series.__getitem__ treating keys as positions is deprecated. In a future version, integer keys will always be treated as labels (consistent with DataFrame behavior). To access a value by position, use `ser.iloc[pos]`
  plot_mr(x_line, wls.params[0] * x_line, merged)
Done