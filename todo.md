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

