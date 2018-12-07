'''
Created on 4 Dec 2018

@author: mate
'''


from cameo import load_model, fba

model = load_model('models/CD4T1670.xml')
new_solution = model.optimize()
fba_result = fba(model)
print(fba_result)
print(fba_result.data_frame)

fba_result.display_on_map("RECON1.Carbohydrate metabolism")