'''
Created on 3 Jul 2019

@author: mate
'''


# import cobra.io

# import os
# 
# model = cobra.test.create_test_model("textbook")
# print(model)
# 
# cobra.io.save_matlab_model(model, "test.mat")
# 
# newModel = cobra.io.load_matlab_model("/home/mate/code/met_model/src/test.mat")
# print(newModel.reactions)

# newModel = cobra.io.load_matlab_model("/home/mate/code/met_model/src/data/WTmodel.mat")# , variable_name="mini_textbook")
# print(newModel)
# # print(newModel.reactions)
# 
# cobra.io.save_json_model(newModel, "WTmodel.json")

from cameo import load_model, fba

model = load_model('WTmodel.json')
print(len(model.genes))
print(len(model.metabolites))
# print(model.reactions)
new_solution = model.optimize()
fba_result = fba(model)
print(fba_result)
print(fba_result.data_frame)
print(len(fba_result.data_frame))
# print(model.reactions)

# for reactionO in model.reactions:
#   for reactionI in reactionO.genes:
#     print(reactionI)

# fba_result.display_on_map("RECON1.Carbohydrate metabolism")

