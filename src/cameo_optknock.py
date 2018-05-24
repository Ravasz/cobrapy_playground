'''
Created on 24 May 2018

@author: mate
'''
from cameo import models

model = models.bigg.iJO1366

wt_solution = model.optimize()
growth = wt_solution.fluxes["BIOMASS_Ec_iJO1366_core_53p95M"]
acetate_production = wt_solution.fluxes["EX_ac_e"]

print(growth, acetate_production)

from cameo.strain_design import OptKnock

optknock = OptKnock(model, fraction_of_optimum=0.1)

result = optknock.run(max_knockouts=1, target="EX_ac_e", biomass="BIOMASS_Ec_iJO1366_core_53p95M")


