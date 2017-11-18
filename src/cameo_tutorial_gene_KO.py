'''
Created on 16 Nov 2017

@author: mate
'''
from cameo import models

model = models.bigg.iJO1366

wt_solution = model.optimize() # model.solve() does not work, but model.optimize() does
growth = wt_solution.fluxes["BIOMASS_Ec_iJO1366_core_53p95M"]
acetate_production = wt_solution.fluxes["EX_ac_e"]

from cameo import phenotypic_phase_plane

p = phenotypic_phase_plane(model, variables=['BIOMASS_Ec_iJO1366_core_53p95M'], objective='EX_ac_e')
p.plot(points=[(growth, acetate_production)])

