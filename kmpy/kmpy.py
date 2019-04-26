import ode_builder

#Reading the filepath and storing them

reactionlist, rateconstantlist, compositionlist = ode_builder.set_paths()


#initializing reactant, product and unique species list
reactants_list = []
products_list = []
species_name = []

for line in open(reactionlist, 'r').readlines():
    reac = ode_builder.Reaction() #Reaction class initiated
    reactants_list.append(reac.getReactantsName(line))
    products_list.append(reac.getProductsName(line))
    current_species = species_name
    #print(current_species)
    species_list = reac.uniqueSpeciesName(line, current_species)
    #print(species_name)
species_list.sort()

print(species_list)
