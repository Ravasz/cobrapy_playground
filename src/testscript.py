'''
Created on 25 May 2017

@author: mate
'''

import cobra

import cobra.test

# "ecoli" and "salmonella" are also valid arguments
model = cobra.test.create_test_model("textbook")
"""
print(len(model.reactions))
print(len(model.metabolites))
print(len(model.genes))

print(model.reactions[29])
"""
for i in model.metabolites:
    print(i)
"""

print(model.reactions.EX_glc__D_e.reversibility)
"""
pgi = model.reactions.get_by_id("ATPS4r") # this is a useful function for handling reactions 

print(pgi.name)
print(pgi.reaction)

pgi.check_mass_balance()

pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
pgi.reaction

atp = model.metabolites.get_by_id("atp_c")
print(atp)

print(atp.name)
print(atp.compartment)
print(model.metabolites.atp_c.formula)

for j in atp.reactions:
    print(j)
    
    
gpr = pgi.gene_reaction_rule
print(gpr)

print(pgi.genes)

for k in pgi.genes:
    print(k)
    
pgi.gene_reaction_rule = "(spam or eggs)"

cobra.manipulation.delete_model_genes(
    model, ["spam"], cumulative_deletions=True)
print("after 1 KO: %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))

cobra.manipulation.delete_model_genes(
    model, ["eggs"], cumulative_deletions=True)
print("after 2 KO:  %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))


pie = model.metabolites.get_by_id("pi_e")
print(pie.name)