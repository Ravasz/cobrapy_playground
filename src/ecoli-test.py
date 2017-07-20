'''
Created on 20 Jul 2017

@author: mate
'''
from cobra.io import load_json_model
import cobra.test

#model = load_json_model("C:/Users/maran/Downloads/iAF1260.json")

from cobra import Model, Reaction, Metabolite

model = cobra.test.create_test_model('textbook')

print(len(model.reactions))
print(len(model.metabolites))
print(len(model.genes))




print(model.reactions)
ex_glu = model.reactions.get_by_id("EX_glc__D_e")
ex_oxygen=model.reactions.get_by_id("EX_o2_e")
print(ex_glu.name)
print(ex_glu.reaction)
ex_glu.reversibility



ex_glu.lower_bound

solution = model.optimize()
print(solution)
print(ex_glu.bounds)
print(ex_oxygen.bounds)
model.summary()

"""
for i in range(0,1000,10):
    with model as model:
        ex_glu.lower_bound = -18.5
        ex_glu.upper_bound = 1
        ex_oxygen.lower_bound = -1*i
        ex_oxygen.upper_bound = 1 
        #print(ex_glu.bounds)
        print(ex_oxygen.bounds)
        model.optimize()
        print(model.objective.value)
        
ex_succ = model.reactions.get_by_id("EX_succ_e")


with model as model:
    ex_succ.lower_bound = -1000
    ex_succ.upper_bound = 1
    ex_glu.lower_bound = 0
    ex_glu.upper_bound = 1
    ex_oxygen.lower_bound = -1000
    ex_oxygen.upper_bound = 1 
    #print(ex_glu.bounds)
    print(ex_oxygen.bounds)
    model.optimize()        
    print(model.objective.value)
    """

aTPM = model.reactions.get_by_id("ATPM")
    
with model as model:
    model.objective = 'ATPM'
    ex_glu.lower_bound = -1
    ex_glu.upper_bound = -1
    ex_oxygen.lower_bound = -1000
    ex_oxygen.upper_bound = 1 
    #print(ex_glu.bounds)
    #print(ex_oxygen.bounds)
    model.optimize()        
    print(model.objective.value)
