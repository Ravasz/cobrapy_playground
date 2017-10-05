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
oxoGlut.subsystem = 'Pyruvate Metabolism'
oxoGlut.lower_bound = 0.  # This is the default
oxoGlut.upper_bound = 1000.  # This is the default
bigg.add_reactions([oxoGlut])

succCoASyn = Reaction("SCoAs")
succCoASyn.name = 'succinyl-CoA synthetase'
succCoASyn.subsystem = 'Pyruvate Metabolism'
succCoASyn.lower_bound = 0.  
succCoASyn.upper_bound = 1000.  
bigg.add_reactions([succCoASyn])

hydroBH = Reaction("4HBD")
hydroBH.name = '4 hydroxybutyrate dehydrogenase'
hydroBH.subsystem = 'Pyruvate Metabolism'
hydroBH.lower_bound = 0.  
hydroBH.upper_bound = 1000.  
bigg.add_reactions([hydroBH])

hydroBCOAT = Reaction("4HBCOAT")
hydroBCOAT.name = '4 hydroxybutyryl-CoA transferase'
hydroBCOAT.subsystem = 'Pyruvate Metabolism'
hydroBCOAT.lower_bound = 0.  
hydroBCOAT.upper_bound = 1000.  
bigg.add_reactions([hydroBCOAT])

hydroBCOAR = Reaction("4HBCOAR")
hydroBCOAR.name = '4 hydroxybutyryl-CoA reductase'
hydroBCOAR.subsystem = 'Pyruvate Metabolism'
hydroBCOAR.lower_bound = 0.  
hydroBCOAR.upper_bound = 1000.  
bigg.add_reactions([hydroBCOAR])


Succ = bigg.metabolites.get_by_id("sucsal_c")
CO2 = bigg.metabolites.get_by_id("co2_c")
AKG = bigg.metabolites.get_by_id("akg_c")
SucCoa = bigg.metabolites.get_by_id("succoa_c")
COA = bigg.metabolites.get_by_id("coa_c")
NADH = bigg.metabolites.get_by_id("nadh_c")
NAD = bigg.metabolites.get_by_id("nad_c")
ACCOA = bigg.metabolites.get_by_id("accoa_c")
AC = bigg.metabolites.get_by_id("ac_c")


HB_c = Metabolite(
    '4HB_c',
    formula='C4H7O3',
    name='4HB_c',
    compartment='c')


HBcoa_c = Metabolite(
    '4HBcoa_c',
    formula='C25H42N7O18P3S',
    name='4HBcoa_c',
    compartment='c')


HBaldehyde = Metabolite(
    'HBaldehyde_c',
    formula='C4H8O2',
    name='HBaldehyde_c',
    compartment='c')


BDO = Metabolite(
    'BDO_c',
    formula='C4H10O2',
    name='BDO_c',
    compartment='c')


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

succCoASyn.add_metabolites({
    SucCoa: -1.0,
    NADH: -1.0,
    COA: 1.0,
    NAD: 1.0,
    Succ: 1.0,
})

hydroBH.add_metabolites({
    Succ: -1.0,
    NADH: -1.0,
    NAD: 1.0,
    HB_c: 1.0,
})

hydroBCOAT.add_metabolites({
    ACCOA: -1.0,
    AC: 1.0,
    HBcoa_c: 1.0,
    HB_c: -1.0,
})

hydroBCOAR.add_metabolites({
    NADH: -1.0,
    NAD: 1.0,
    COA: 1.0,
    HB_c: -1.0,
    HBaldehyde: 1.0,
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

