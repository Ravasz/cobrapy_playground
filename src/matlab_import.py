'''
Created on 3 Jul 2019

@author: mate
'''


# import cobra.io
# #   
# # import os
# #  
# # model = cobra.test.create_test_model("textbook")
# # print(model)
# #  
# # cobra.io.save_matlab_model(model, "test.mat")
#  
# # newModel = cobra.io.load_matlab_model("/home/mate/code/met_model/src/test.mat")
# #print(newModel.reactions)
# print("test")
# newModel = cobra.io.load_matlab_model("/home/mate/code/met_model/src/data/KOmodel.mat")# , variable_name="mini_textbook")
# print(newModel)
# # print(newModel.reactions)
#  
# cobra.io.save_json_model(newModel, "KOmodel.json")

def write_results():
  """optimize results and write out to a file"""
  from cameo import load_model, fba
  
  outF = open("KOFBAFluxes-2.txt","w")
  
  model = load_model('KOmodel.json')
  print(len(model.genes))
  print(len(model.metabolites))
  print(len(model.reactions))
  # new_solution = model.optimize()
  fba_result = fba(model)
  print(fba_result)
  # print(fba_result.data_frame)
  print(len(fba_result.data_frame))
  # print(model.reactions)
    
  print(dir(fba_result))
  print(fba_result.data_frame)
  fbaD = fba_result.data_frame.T.to_dict()
  fbaCorrD = {}
  for keyS, valueF in fbaD.items():
    fbaCorrD[keyS] = valueF["flux"]
 
  for keyS, valueF in fbaCorrD.items():
    # print(keyS, valueF)
    outF.write(str(keyS) + ", " + str(valueF) + "\n")
    # print(reactionI)
#       
  
    
  # fba_result.display_on_map("RECON1.Carbohydrate metabolism")
  

def find_gene():
  """create optimized model from json file and find genes and reactions in it"""
  from cameo import load_model
  
  # [('HGNC:1093', 0.8281729165648598), ('HGNC:8889', 0.8281729165648598), ('HGNC:8888', 0.8281729165648598), ('HGNC:8898', 0.8281729165648598), ('HGNC:8896', 0.8281729165648598), 
  # ('HGNC:16270', 0.6949058320029933), ('HGNC:8124', 0.6), ('HGNC:3700', 0.6), ('HGNC:4336', 0.6), ('HGNC:4335', 0.6)]
  
  model = load_model('WTmodel.json')
  print(len(model.genes))
  print(len(model.metabolites))
  print(len(model.reactions))
  
  model.optimize()
  
  for modelO in model.reactions:
    # print(dir(modelO))
    # print(modelO.name, modelO.reaction, modelO.nice_id, modelO.annotation)
    if str(modelO.nice_id) == "PGK":
      print(modelO.name, modelO.reaction, modelO.nice_id, modelO.annotation, modelO.flux)
    elif str(modelO.nice_id) == "O2t":
      print(modelO.name, modelO.reaction, modelO.nice_id, modelO.annotation, modelO.flux)
    

    # r1144
# WT
# o2 transport (diffusion) o2[e] <=> o2[c] O2t {} 1.3990599985415846
# Major Facilitator(MFS) TCDB:2.A.18.6.3 glu_L[e] + na1[e] --> glu_L[c] + na1[c] r1144 {} 0.6929768192572117


# KO
# o2 transport (diffusion) o2[e] <=> o2[c] O2t {} 0.7490861815722542
# Major Facilitator(MFS) TCDB:2.A.18.6.3 glu_L[e] + na1[e] --> glu_L[c] + na1[c] r1144 {} 0.0
#



def compare_fluxes():
    
  '''
  Created on 8 Dec 2019
  
  @author: mate
  
  compare WT and KO fluxes made by the script above
  '''
  
  
  KOFile = open("KOFBAFluxes.txt","r")
  
  KODict = {}
  
  for KOLine in KOFile:
    # print(KOLine)
    KOList = KOLine.strip("\n").split(",")
    KODict[KOList[0]] = float(KOList[1])
  
  KOFile.close()
  
  WTFile = open("WTFBAFluxes.txt","r")
  
  WTDict = {}
  
  for WTLine in WTFile:
    WTList = WTLine.strip("\n").split(",")
    WTDict[WTList[0]] = float(WTList[1])
    
  WTFile.close()
  
  mergeD = {}
  
  for idI in KODict:
    mergeD[idI] = WTDict[idI] - KODict[idI]
    
  sortedL = sorted(mergeD.items(), reverse = True, key=lambda kv: kv[1])
  
  # print(sortedL[:10])
  rowCount = 0
  for sortI in sortedL:
    print(sortI[0], sortI[1])
    rowCount += 1
    if rowCount == 20: break


# write_results()
# compare_fluxes()
find_gene()
  