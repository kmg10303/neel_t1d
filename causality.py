from dowhy import CausalModel


model = CausalModel(
    data=df,
    treatment="IL6"/"beta_exp",
    outcome="T1D"/"bet_out",
    graph="digraph { IL6 -> T1D; G -> T1D; U -> IL6; U -> T1D; }",
    instruments=["Z"]
)

# Refutation tests:
# Add a random common cause. (weight, exercise, diet)
# - refute_estimate("add_unobserved_common_cause")
# 
# Placebo (there's no connection between two variables):
# - refute_estimate("placebo_treatment_refuter")
# 
# 
# 
