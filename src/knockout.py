'''
Created on 13 Sep 2017

@author: mate
'''
import cobra.test

"""
outF = open("doubleKO_PGI.csv","w")
outF.write("Reaction1,Reaction2,objective\n")
"""
model = cobra.test.create_test_model('textbook')
model.optimize()
print(model.objective.value)

mylist = model.reactions
newList = model.reactions



print(len(mylist))
lineCount = 0
rowCount = 0

newObjective = model.reactions.get_by_id("ACALD")

exGlu = model.reactions.get_by_id("EX_glc__D_e")

exGlu.lower_bound = -1000
exGlu.upper_bound = -1000

print(newObjective)

for reaction in mylist:
  print(reaction)
  
  for newReaction in newList:
    if reaction == newReaction: 
      #print("two reactions are identical")
      continue
    #print(newReaction)
    lineCount += 1
    with model as newModel:
      model.objective = newObjective
      """
      if lineCount == 100:
        lineCount = 0
        print(".", end = "")
        rowCount += 1
        if rowCount == 100:
          rowCount = 0
          print()
      """
      reaction.knock_out()
      newReaction.knock_out()
      newModel.optimize()
      #outF.write(str(reaction) + "," + str(newReaction) + "," + str(newModel.objective.value) + "\n")
      if round(float(newModel.objective.value),6) != 0.873922:
        print('%s,%s blocked, new growth rate %f' %(reaction.id, newReaction.id, newModel.objective.value))
          
        

