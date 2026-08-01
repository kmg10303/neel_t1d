To do for 8/8:
- Methodology is slightly different now. This is the final version. Can do a quick write-up on the methodology this week.
- Using GAI, attempt to add the real estimate section of the causal graph. 


To do for 5/31:
- Methodology up to this point (keep 1 paragraph)
- Editing intro and data source section
- python dowhy --> reading through the intro documentation


Implementation Plan:

Task 1: Add more models. 
- WLS w/ different intercepts.
- Median, weighted median (medians of our Wald ratios)

Task 2: Fix the scatter plot.

Task 3: Turn our models into causal models.
- use python dowhy package





**For Neel to research over the next few weeks**

Write-up:
- How have past researchers shown that HLA-DQ1 is a factor for risk of T1D?
- Is HLA-DQ1 this a good candidate for an MR study?
    - Too complex, lots of confounding results.
    - Protein could have a lot of different exposures. 
        - rs9272346
- Data sources we can write out. Look into the Sun BB IL2B study. 
- Expected outcomes. Do we already know that IL2B might increase risk for t1D? If we do, we probably expect our result to show this. 
    - As evidence for your expectations, you should use the DAG.

- Create a final hypothesis. ("Increase in IL2B proteins increases risk for T1D"... "Increase in IL2B proteins decrease/increases X/Y which, in turn, increases risk for T1D").

How to run function:

- get API key:
    - https://api.opengwas.io/
    - Sign in with github or microsoft
    - Create and copy your token
    - Create file called .ieugwaspy.json
    - In that file paste in this:
        - {"base_url": "https://api.opengwas.io/api/", "jwt": [INSERT TOKEN]}

. venv/bin/activate
python mr_t1d.py

Code:
- Add more SNPs from the dataset we found in mr_t1d.py

