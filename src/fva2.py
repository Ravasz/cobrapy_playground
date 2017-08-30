'''
Created on 30 Aug 2017

@author: mate

This is a test for FVA. 
This can show how FVA can give a range of solutions based on preset constraints
'''

import cobra.test

model = cobra.test.create_test_model("textbook")

from cobra.flux_analysis import flux_variability_analysis

print (flux_variability_analysis(model, model.reactions))

# model.objective = "ME1"

ex_glu = model.reactions.get_by_id("EX_glc__D_e")
ex_glu.upper_bound = 1
ex_glu.lower_bound = 0

succinate = model.reactions.get_by_id("EX_succ_e")
succinate.upper_bound = 1
succinate.lower_bound = -20

print (flux_variability_analysis(model, model.reactions))

ex_O = model.reactions.get_by_id("EX_o2_e")
ex_O.upper_bound = -17
ex_O.lower_bound = -17

succinate.lower_bound = 0



mylist = []

for i in range(20):
  ex_glu.upper_bound = i * (-1)
  ex_glu.lower_bound = i * (-1)
  solu = model.optimize()
  mylist.append(solu.objective_value)
  # print(solu.objective_value)
  print(solu.reactions.get_by_id("ACALD").flux)
  
print(mylist)  


