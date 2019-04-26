import ode_builder

#Reading the filepath and storing them
reactionlist, rateconstantlist, compositionlist = ode_builder.set_paths()

#builiding the reactant, product and unique species list
reactant_list, product_list, unique_species = \
ode_builder.build_species_list(reactionlist)

print(unique_species)
