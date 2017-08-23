'''
Created on 27 Jul 2017

@author: mate
'''

"""
from escher import Builder
b = Builder(map_name="iJO1366.Central metabolism")
b.display_in_browser()
"""

"""
import cobra, cobra.test

from escher import Builder

newModel = cobra.test.create_test_model('textbook')

with newModel as newModel:
    newModel.optimize()
    
cobra.io.write_sbml_model(newModel, "textbook_export.xml")
my_cobra_model = cobra.io.read_sbml_model('textbook_export.xml')
b = Builder(model=my_cobra_model)
b.display_in_browser()
"""
import escher
import escher.urls
import cobra
import cobra.test
import json
import os

"""
d = escher.urls.root_directory
print('Escher directory: %s' % d)

escher.list_available_maps()
"""
"""
b = escher.Builder(map_name='e_coli_core.Core metabolism')
#b.display_in_browser()

model = cobra.io.json.from_json(escher.plots.model_json_for_name('e_coli_core'))
sol = model.optimize()
sol.fluxes.to_json("solution2.json")
jsonFile = open("solution2.json","r")
jsonItem = jsonFile.read()
testDict = json.loads(jsonItem)
print(testDict)

print('Growth rate: %.2f' % sol.f)
"""

import cobra.test
model = cobra.test.create_test_model("textbook")
ex_oxygen = model.reactions.get_by_id("EX_o2_e")
with model as newModel:
    ex_oxygen.lower_bound = 0
    ex_oxygen.upper_bound = 1
    newSolution = newModel.optimize()
    
    newSolution.fluxes.to_json("solution3.json")
    jsonFile = open("solution3.json","r")
    jsonItem = jsonFile.read()
    testDict = json.loads(jsonItem)
    print(testDict)
    
    

    
    
    b = escher.Builder(map_name='e_coli_core.Core metabolism',
                       reaction_data=testDict,
                       # change the default colors
                       reaction_scale=[{'type': 'min', 'color': '#cccccc', 'size': 4},
                                       {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                       {'type': 'max', 'color': '#ff0000', 'size': 40}],
                       # only show the primary metabolites
                       hide_secondary_metabolites=False)
    #b.display_in_browser(js_source='web')
    b.display_in_browser()