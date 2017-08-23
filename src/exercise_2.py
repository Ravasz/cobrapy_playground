'''
Created on 23 Aug 2017

@author: mate
'''

import cobra.test
from cobra import Reaction

testModel = cobra.test.create_test_model('textbook')

with testModel as currModel:
  currModel.objective = "ATPM"
  ex_glu = testModel.reactions.get_by_id("EX_glc__D_e")
  ex_glu.lower_bound = -1
  ex_glu.upper_bound = -1
  
  testModel.optimize()
  print(testModel.objective.value)
  
  nadReaction = Reaction('nadhd')
  nadReaction.name = 'nadh drain'
  nadReaction.lower_bound = 0.  # This is the default
  nadReaction.upper_bound = 1000.  # This is the default
  nadh_c = currModel.metabolites.get_by_id("nadh_c")
  nad_c = currModel.metabolites.get_by_id("nad_c")
  h_c = currModel.metabolites.get_by_id("h_c")
  
  nadReaction.add_metabolites({
    nadh_c: -1.0,
    nad_c: 1.0,
    h_c: 1.0,
    })

  currModel.add_reactions([nadReaction])
  
  # print(currModel.reactions)
  nadhd = currModel.reactions.get_by_id("nadhd")
  # print (nadhd.metabolites)
  
  currModel.objective = "nadhd"
  currModel.optimize()
  print(currModel.objective.value)