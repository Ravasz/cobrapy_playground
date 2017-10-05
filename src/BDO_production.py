'''
Created on 5 Oct 2017

@author: mate
'''
from cobra import Model, Reaction, Metabolite

from cobra.io import load_json_model

bigg = load_json_model("/home/mate/code/cobrapy-playground/src/models/iJR904.json")

print(len(bigg.reactions))
print(len(bigg.metabolites))
print(len(bigg.genes))

# print(bigg.metabolites)

oxoGlut = Reaction('2OxDe')
oxoGlut.name = '2 oxoglutarate decarboxylase'
oxoGlut.subsystem = ''
oxoGlut.lower_bound = 0.  # This is the default
oxoGlut.upper_bound = 1000.  # This is the default
bigg.add_reactions([oxoGlut])

succCoASyn = Reaction("SCoAs")
succCoASyn.name = 'succinyl-CoA synthetase'
succCoASyn.subsystem = ''
succCoASyn.lower_bound = 0.  
succCoASyn.upper_bound = 1000.  
bigg.add_reactions([succCoASyn])



Succ = bigg.metabolites.get_by_id("sucsal_c")
CO2 = bigg.metabolites.get_by_id("co2_c")
AKG = bigg.metabolites.get_by_id("akg_c")
SucCoa = bigg.metabolites.get_by_id("succoa_c")
COA = bigg.metabolites.get_by_id("coa_c")
NADH = bigg.metabolites.get_by_id("nadh_c")
NAD = bigg.metabolites.get_by_id("nad_c")


oxoGlut.add_metabolites({
    Succ: 1.0,
    CO2: 1.0,
    AKG: -1.0,
})

succCoASyn.add_metabolites({
    SucCoa: -1.0,
    NADH: -1.0,
    COA: 1.0,
    NAD: 1.0,
    Succ: 1.0,
})



"""
Ssem_c = Metabolite(
    'Ssem_c',
    formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein',
    compartment='c')
co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
malACP_c = Metabolite(
    'malACP_c',
    formula='C14H22N2O10PRS',
    name='Malonyl-acyl-carrier-protein',
    compartment='c')
h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
ddcaACP_c = Metabolite(
    'ddcaACP_c',
    formula='C23H43N2O8PRS',
    name='Dodecanoyl-ACP-n-C120ACP',
    compartment='c')
"""


print(oxoGlut.metabolites)


4HB_c = Metabolite(
    '4HB_c',
    formula='C4H7O3',
    name='4HB_c',
    compartment='c')


4HBcoa_c = Metabolite(
    '4HBcoa_c',
    formula='C25H42N7O18P3S',
    name='4HBcoa_c',
    compartment='c')
