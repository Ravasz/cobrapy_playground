'''
Created on 15 Jun 2017

@author: mate
'''

import cobra.test
model = cobra.test.create_test_model("textbook")

solution = model.optimize()

print(solution.objective_value)
# print(solution.shadow_prices)

print(model.summary())

model.metabolites.nadh_c.summary()

model.metabolites.atp_c.summary()

# atpmO = model.reactions.get_by_id("ATPM")
# growthO = model.reactions.get_by_id("Biomass_Ecoli_core")
# pfkO = model.reactions.get_by_id("PFK")

# print("\nATPM flux when optimized for growth: " + str(atpmO.flux))
# print("Growth when optimized for growth: " + str(growthO.flux))
# print("PFK flux when optimized for growth: " + str(pfkO.flux))


from cobra.util.solver import linear_reaction_coefficients
print(linear_reaction_coefficients(model))

# change the objective to ATPM
model.objective = "ATPM"

# The upper bound should be 1000, so that we get
# the actual optimal value
model.reactions.get_by_id("ATPM").upper_bound = 1000.
linear_reaction_coefficients(model)
model.optimize().objective_value

# print("\nATPM flux when optimized for ATPM: " + str(atpmO.flux))
# print("Growth when optimized for ATPM: " + str(growthO.flux))
# print("PFK flux when optimized for ATPM: " + str(pfkO.flux))




# model.objective = "PFK"

# The upper bound should be 1000, so that we get
# the actual optimal value
# model.reactions.get_by_id("PFK").upper_bound = 1000.
# linear_reaction_coefficients(model)
# model.optimize()

# print("\nATPM flux when optimized for PFK: " + str(atpmO.flux))
# print("Growth when optimized for PFK: " + str(growthO.flux))
# print("PFK flux when optimized for PFK: " + str(pfkO.flux))

from cobra.flux_analysis import flux_variability_analysis

print(flux_variability_analysis(model, model.reactions[:10]))

print(cobra.flux_analysis.flux_variability_analysis(model, model.reactions[:10], fraction_of_optimum=0.9))

print("-------")

loop_reactions = [model.reactions.FRD7, model.reactions.SUCDi]
print(flux_variability_analysis(model, reaction_list=loop_reactions, loopless=True))

loop_reactions = [model.reactions.FRD7, model.reactions.SUCDi]
print(flux_variability_analysis(model, reaction_list=loop_reactions, loopless=False))
