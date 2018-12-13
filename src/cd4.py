'''
Created on 4 Dec 2018

@author: mate

load in CD4 model, extract all proteins and save them to a new .csv
'''


from cameo import load_model, fba

model = load_model('models/CD4T1670.xml')
new_solution = model.optimize()
fba_result = fba(model)
print(fba_result)
print(fba_result.data_frame)
print(len(fba_result.data_frame))
# print(model.reactions)

for reactionO in model.reactions:
  for reactionI in reactionO.genes:
    print(reactionI)

fba_result.display_on_map("RECON1.Carbohydrate metabolism")

