"""
Created on 1 Jun 2017

@author: mate

import cobra, sys

print(cobra.__version__)

print(sys.version)
"""
import cobra.test



model = cobra.test.create_test_model('textbook')
print(model.objective.value)
print(model.reactions[49])
GLT = model.reactions.get_by_id("ATPM")
print(GLT)
print(GLT.name)
print(GLT.reaction)
mylist = model.reactions[:27] + model.reactions[28:49] + model.reactions[50:]
print(len(mylist))
for reaction in mylist:
    with model as model:
        reaction.knock_out()
        model.optimize()
        if model.objective.value > 0:
            print('%s blocked (bounds: %s), new growth rate %f' %
                  (reaction.id, str(reaction.bounds), model.objective.value))

"""
print('original objective: ', model.objective.expression)
with model:
    model.objective = 'ATPM'
    print('print objective in first context:', model.objective.expression)
    with model:
        model.objective = 'ACALD'
        print('print objective in second context:', model.objective.expression)
    print('objective after exiting second context:',
          model.objective.expression)
print('back to original objective:', model.objective.expression)
            
print("Reactions")
print("---------")
for x in model.reactions:
    print("%s : %s" % (x.id, x.reaction))

print("")
print("Metabolites")
print("-----------")
for x in model.metabolites:
    print('%9s : %s' % (x.id, x.formula))

print("")
print("Genes")
print("-----")
for x in model.genes:
    associated_ids = (i.id for i in x.reactions)
    print("%s is associated with reactions: %s" %
          (x.id, "{" + ", ".join(associated_ids) + "}"))                
            
"""