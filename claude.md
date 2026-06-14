Project Summary: Causal AI × T1D MR Analysis
Core question: Do genetically predicted cytokine levels (specifically IL-6, with IL-1β as a comparison) causally influence Type 1 Diabetes risk?
What exists:

A working ieugwaspy + OpenGWAS pipeline that queries SNP associations, harmonizes instruments, and runs IVW regression via WLS
Multi-SNP analysis (3–4 instruments) producing a causal estimate, F-statistics, R², and p-values
A working matplotlib scatter plot showing the regression line (though plt.errorbar() for SNP data points was still missing at last session)
A draft writeup that mentions "causal AI" but hasn't yet concretely implemented any causal AI components

Key known issues:

tophits returns empty — most likely an expired JWT token on the OpenGWAS API (expired around March 3)
The best fix is switching exposure GWAS to Ahola-Olli et al. 2017 (ieu-b- batch, 41 cytokines, N=8,293) or Zhao et al. (91-cytokine pQTL, GCST90274758–GCST90274848, N=14,824)
IL-6 is the strongest primary exposure — rs4537545 is a very well-powered instrument
The IL-1β run (beta=−27.7, p=0.243, F=1.87) is a weak instrument result and should be kept as a contrast, not discarded
The pipeline still lacks MR-Egger, weighted median, and sensitivity analyses
The "causal AI" framing is currently a hand-wave — no causal discovery algorithms have been implemented yet

