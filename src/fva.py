'''
Created on 22 Jun 2017

@author: mate
'''

import cobra.test
model = cobra.test.create_test_model("textbook")

solution = model.optimize()

"""
print(solution)

solution.objective_value

model.summary()

model.metabolites.nadh_c.summary()

model.metabolites.atp_c.summary()

biomass_rxn = model.reactions.get_by_id("Biomass_Ecoli_core")

"""

from cobra.util.solver import linear_reaction_coefficients
linear_reaction_coefficients(model)

# change the objective to ATPM
model.objective = "ATPM"

# The upper bound should be 1000, so that we get
# the actual optimal value
model.reactions.get_by_id("ATPM").upper_bound = 1000.
linear_reaction_coefficients(model)

model.optimize().objective_value

from cobra.flux_analysis import flux_variability_analysis

print(flux_variability_analysis(model, model.reactions))

cobra.flux_analysis.flux_variability_analysis(
    model, model.reactions[:10], fraction_of_optimum=0.9)



"""
loop_reactions = [model.reactions.FRD7, model.reactions.SUCDi]
flux_variability_analysis(model, reaction_list=loop_reactions, loopless=False)

flux_variability_analysis(model, reaction_list=loop_reactions, loopless=True)

model.optimize()
model.summary(fva=0.95)

model.metabolites.pyr_c.summary(fva=0.95)


model.objective = 'Biomass_Ecoli_core'
fba_solution = model.optimize()
pfba_solution = cobra.flux_analysis.pfba(model)
print(pfba_solution)

print(abs(fba_solution.fluxes["Biomass_Ecoli_core"] - pfba_solution.fluxes[
    "Biomass_Ecoli_core"]))

"""