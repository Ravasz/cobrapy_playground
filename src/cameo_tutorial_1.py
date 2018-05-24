'''
Created on 8 Nov 2017

@author: mate
'''

print("Hello Cameo")

from cameo import load_model
model = load_model('iJO1366') 

from cameo import fba
fba_result = fba(model)
print(fba_result)
print(fba_result.data_frame)
 
# fba_result.display_on_map("iJO1366.Central metabolism")

from cameo import pfba
pfba_result = pfba(model)
print(pfba_result.data_frame)
print(pfba_result.objective_value)
print(pfba_result[model.reactions.BIOMASS_Ec_iJO1366_core_53p95M])

print(model.reactions.PGI.upper_bound)
print(model.reactions.PGI.lower_bound)

model.reactions.PGI.knock_out()

pfba_knockout_result = pfba(model)

print(model.reactions.PGI.upper_bound)
print(model.reactions.PGI.lower_bound)

print(pfba_knockout_result.objective_value)
print(pfba_knockout_result[model.reactions.BIOMASS_Ec_iJO1366_core_53p95M])

from cameo.strain_design import OptKnock




# from cobra import io

# io.save_json_model(model,"test.json")



#from cameo import models
#biggModelList = models.index_models_bigg()
#newModel = models.bigg.iJN746
#print("model imported")
#print(newModel.reactions)

#from cameo import api
#api.design(product="vanillin")

# from optlang import Model, Variable, Constraint, Objective
# 
# # All the (symbolic) variables are declared, with a name and optionally a lower and/or upper bound.
# x1 = Variable('x1', lb=0)
# x2 = Variable('x2', lb=0)
# x3 = Variable('x3', lb=0)
# 
# # A constraint is constructed from an expression of variables and a lower and/or upper bound (lb and ub).
# c1 = Constraint(x1 + x2 + x3, ub=100)
# c2 = Constraint(10 * x1 + 4 * x2 + 5 * x3, ub=600)
# c3 = Constraint(2 * x1 + 2 * x2 + 6 * x3, ub=300)
# 
# # An objective can be formulated
# obj = Objective(10 * x1 + 6 * x2 + 4 * x3, direction='max')
# 
# # Variables, constraints and objective are combined in a Model object, which can subsequently be optimized.
# model = Model(name='Simple model')
# model.objective = obj
# model.add([c1, c2, c3])
# 
# status = model.optimize()
# 
# print("status:", model.status)
# print("objective value:", model.objective.value)
# print("----------")
# for var_name, var in model.variables.iteritems():
#     print(var_name, "=", var.primal)