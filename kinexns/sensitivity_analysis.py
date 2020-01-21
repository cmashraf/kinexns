from rdkit import Chem
import rdkit.Chem.Draw as Draw
from PIL import Image, ImageOps
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def gen_params(sample_number, sa_path, input_file_name, output_file_name):
    """
    Generate the parameter sets for the Sobol sensitivity analysis.
    Saves a file with the information required for the analysis
    that will be performed later.
    Parameters
    -------------
    sample_number   : number of samples to be generated. Total of parameters
                    generated will be N * (2D+2) where N = sample number and
                    D = number of parameters to vary which is listed in
                    input_file_name
    sa_path         : path of the files to read and write from
    input_file_name : name of the input file located at sa_path. This file
                    should have three columns - parameter name, lower limit
                    and upper limit
    output_file_name    : name at the output file where the parameter
                        combiantions should be writtten located at
                        sa_path
    Returns
    --------
    param_sets      : numpy ndarray written at output_file_name
                    an ndarray where each row is one set of parameter
                    values.  You must run your model (in whatever environment
                    is appropriate) with each of the sets of parameters from
                    this array.  The output must be stored in the same order
                    as given in this parameter set array (one row of results
                    for each row of parameters).
    """

    input_file = sa_path + input_file_name
    output_file = sa_path + output_file_name
    os.system('python -m SALib.sample.saltelli -n {} -p {} -o '
              '{}'.format(sample_number, input_file, output_file))


def analyze_sensitivity(path, column, delimiter, order, name,
                        parallel=False, processors=4):
    """
    Perform the sensitivity analysis after you have run your model
    with all the parameters from gen_params().  This is done from
    the command line because it is faster and gives the option to
    specify the column of the results file to analyze.  Parallel
    processing is possible.  Results are saved to a file using the
    name parameter.
    Parameters
    ----------
    path        : str
                 the path to the saparams* file that contains
                 the problem definition.
    column     : int
                 integer specifying the column number of the results to
                 analyze (zero indexed).
    delimiter  : str
                 string specifying the column delimiter used in the results.
    order      : int
                 the maximum order of sensitivity indices [1 or 2].
    name       : str
                 the name of the output measure to use when saving
                 the sensitivity analysis results to a file.
    parallel   : bool, optional
                 boolean indicating whether to use parallel processing.
    processors : int, optional
                 if parallel is True, this is an integer specifying the number
                 of processors to use.
    Returns
    --------
    A file named as 'analysis_{}.txt with sensitivity analysis results
    for that output measure.
    """
    out_file = path + 'analysis_{}.txt' .format(name)
    params_file = path + 'params.txt'
    result = path + 'model_solutions.txt'
    if parallel:
        os.system('python -m SALib.analyze.sobol -p {} -Y {} -c {} '
                  '--delimiter {} --max-order {} --parallel --processors'
                  ' {} > {}' .format(params_file, result, column, delimiter,
                                     order, processors, out_file))
    else:
        os.system('python -m SALib.analyze.sobol -p {} -Y {} -c {} '
                  '--delimiter {} --max-order {} > {}'
                  .format(params_file, result, column, delimiter,
                          order, out_file))


def default_draw_options():
    """
   This function returns an RDKit drawing options object with
   default drawing options. No argument needed
    """

    opts = Draw.DrawingOptions()
    # opts.elemDict = defaultdict(lambda: (0,0,0)) # all atoms are black
    opts.noCarbonSymbols = True
    opts.selectColor = (1, 0, 0)
    opts.wedgeBonds = True
    return opts


def mol_to_image(mol, max_size=(1000, 1000), kekulize=True, options=None,
                 canvas=None, **kwargs):
    """
    Wrapper for RDKit's MolToImage. If mol == None, an arrow is drawn
    Parameters
    -----------
    mol     : str
            the smiles string of the molecule to be drawn.
            If mol == None, an arrow is drawn.
    max_size    : tuple
                maximum size in pixel of the image to be drawn
                default (1000, 1000)
    kekulize    : boolean
                run kekulization routine on input mol (default True)
    options     : None
    canvas      : None
    Returns
    --------------
    A rdkit image of the molecule
    """

    if not options:
        options = default_draw_options()
    if mol == '->':
        sub_img_size = (160, 160)
        img, canvas = Draw._createCanvas(sub_img_size)
        p0 = (10, sub_img_size[1] // 2)
        p1 = (sub_img_size[0] - 10, sub_img_size[1] // 2)
        p3 = (sub_img_size[0] - 20, sub_img_size[1] // 2 - 10)
        p4 = (sub_img_size[0] - 20, sub_img_size[1] // 2 + 10)
        print('{}'.format(p1))
        canvas.addCanvasLine(p0, p1, lineWidth=2, color=(0, 0, 0))
        canvas.addCanvasLine(p3, p1, lineWidth=2, color=(0, 0, 0))
        canvas.addCanvasLine(p4, p1, lineWidth=2, color=(0, 0, 0))
        if hasattr(canvas, 'flush'):
            canvas.flush()
        else:
            canvas.save()
        return img
    elif mol is not None:
        return Draw.MolToImage(mol, size=max_size, kekulize=kekulize,
                               options=options, canvas=canvas, **kwargs)
    else:  # retro arrow or error
        sub_img_size = (80, 80)
        (a, b) = sub_img_size
        img, canvas = Draw._createCanvas(sub_img_size)
        canvas.addCanvasLine((10, b // 2 - 7), (a - 17, b // 2 - 7),
                             lineWidth=2, color=(0, 0, 255))
        canvas.addCanvasLine((10, b // 2 + 7), (a - 17, b // 2 + 7),
                             lineWidth=2, color=(0, 0, 255))
        canvas.addCanvasLine((a - 24, b // 2 - 14), (a - 10, b // 2),
                             lineWidth=2, color=(0, 0, 255))
        canvas.addCanvasLine((a - 24, b // 2 + 14), (a - 10, b // 2),
                             lineWidth=2, color=(0, 0, 255))
        if hasattr(canvas, 'flush'):
            canvas.flush()
        else:
            canvas.save()
        return img


def trim_img_by_white(img, padding=0):
    """
    This function takes a PIL image, img, and crops it to the minimum
    rectanglebased on its whiteness/transparency.
    5 pixel padding used automatically.
    """

#   Convert to array
    as_array = np.array(img)  # N x N x (r,g,b,a)
#    print(as_array)

# Set previously-transparent pixels to white
    # as_array[as_array[:, :, 3] == 0] = [255, 255, 255]

# Content defined as non-white and non-transparent pixel
    has_content = np.sum(as_array, axis=2, dtype=np.uint32) != 255 * 3
    xs, ys = np.nonzero(has_content)

# Crop down
    x_range = max([min(xs) - 5, 0]), min([max(xs) + 5, as_array.shape[0]])
    y_range = max([min(ys) - 5, 0]), min([max(ys) + 5, as_array.shape[1]])

    as_array_cropped = as_array[x_range[0]:x_range[1],
                                y_range[0]:y_range[1], 0:3]

    img = Image.fromarray(as_array_cropped, mode='RGB')

    return ImageOps.expand(img, border=padding, fill=(255, 255, 255, 0))


def stitch_pils_horizontally(imgs):
    """
    This function takes a list of PIL images and concatenates
    them onto a new image horizontally, with each one
    vertically centered.
    """


# Create blank image (def: transparent white)
    heights = [img.size[1] for img in imgs]
    height = max(heights)
    widths = [img.size[0] for img in imgs]
    width = sum(widths)

    res = Image.new('RGBA', (width, height), (255, 255, 255, 255))

# Add in sub-images
    for i, img in enumerate(imgs):
        offset_x = sum(widths[:i])  # left to right
        offset_y = (height - heights[i]) / 2
        res.paste(img, (int(offset_x), int(offset_y)))

    return res


def stitch_pils_vertically(imgs):
    """
    This function takes a list of PIL images and concatenates
    them onto a new image vertically, with each one
    horizontally centered.
    """

# Create blank image (def: transparent white)
    heights = [img.size[1] for img in imgs]
    height = sum(heights)
    widths = [img.size[0] for img in imgs]
    width = max(widths)

    res = Image.new('RGBA', (width, height), (255, 255, 255, 255))

# Add in sub-images
    for i, img in enumerate(imgs):
        offset_x = (width - widths[i]) / 2  # left to right
        offset_y = sum(heights[:i])
        res.paste(img, (int(offset_x), int(offset_y)))

    return res


def reaction_to_image(rxn, kekulize=True,
                      options=None, **kwargs):
    """
    Modification of RDKit's ReactionToImage to allow for each molecule
    to have a different drawn size. rxn is an RDKit reaction object
    warning: this function adds hydrogens as it sees fit
    Parameters
    ------------
    rxn     : str
             A raction string generated from build reaction function
    kekulize   : boolean
                run kekulization routine on input mol (default True)
    options     : None
    Returns
    --------
    rdkit images of the reactions generated using three separate functions
    trim_img_by_white, mol_to_image and stitch_pils_horizontally
    """

    # Extract mols from reaction
    mols = []
    for item in rxn:
        if item == '':
            mols.append(None)
        else:
            mol = Chem.MolFromSmiles(item)
            mols.append(mol)

#    Generate images for all molecules/arrow
    imgs = [trim_img_by_white(mol_to_image(
            mol, kekulize=kekulize, options=options), padding=5)
            for mol in mols]

# Combine
    return stitch_pils_horizontally(imgs)


def parse_sa_analysis(path):
    """"
    This function reads and processes all the sensitivity analysis results
    in a specified folder and returns a dictionary with the corresponding
    dataframes for first/total order sensitivity indices and second order
    indices (if present).
    Sensitivity analysis results should be in the default SALib output
    format and must start with the word 'analysis'.
    Parameters
    -----------
    path : str, optional
           String containing the relative or absolute path of the directory
           where analysis_*.txt files are stored.  There cannot be any
           files or folders within this directory that start with 'analysis'
           except those generated by the SALib sensitivity analysis.  All
           `analysis*` files in this path should correspond to outputs from
           one sensitivity analysis project, and if second order sensitivity
           indices are included in any of the files they should be present in
           all the others.
    Returns
    --------
    sens_dfs : dict
               Dictionary where keys are the names of the various output
               measures (one output measure per analysis file in the folder
               specified by path).  Dictionary values are a list of pandas
               dataframes.
               sens_dfs['key'][0] is a dataframe with the first and total
               order indices of all the parameters with respect to the "key"
               output variable.
               sens_dfs['key'][1] is a dataframe with the second order
               indices for pairs of parameters (if second order indices are
               present in the analysis file).  If there are no second order
               results in the analysis file then this value is a boolean,
               False.
    """

    filenames = [filename for filename in os.listdir(
                     path) if filename.startswith('analysis')]
    print(filenames)

    dict_dfs = {}

    for filename in filenames:
        name = filename[9:].replace('.txt', '')
# print(name)

        with open(path + filename) as result:
            contents = list([])
            contents.append(result.readlines())
#        print(contents)
# find the line number in the file where 2nd order results appear
            for j, line in enumerate(contents[0]):

                # End this loop when you reach the line that separates
                # the first/total indices from the second order indices
                if line.startswith('\n'):
                    break
# If no second order indices in file
                else:
                    j = False


# If there are second order indices in the file
            if j:
                dict_dfs[name] = [pd.read_csv(path + filename, sep=' ',
                                  nrows=(j - 1)),
                                  pd.read_csv(path + filename, sep=' ',
                                  skiprows=j)]
            else:
                dict_dfs[name] = [pd.read_csv(path + filename, sep=' '),
                                  False]
    return dict_dfs


def get_top_ones(path, species, coeff='total', number=5):
    """
    This function returns the top five reactions with highest absolute
    sensitivity index values.
    Parameters
    -----------
    path    : str, optional
            String containing the relative or absolute path of the directory
            where analysis_*.txt files are stored.  There cannot be any
            files or folders within this directory that start with 'analysis'
            except those generated by the SALib sensitivity analysis.  All
            analysis*` files in this path should correspond to outputs from
            one sensitivity analysis project, and if second order sensitivity
            indices are included in any of the files they should be present in
            all the others.
    species : str
            name of the species/measurable output on whom the analysis will
            be performed
    coeff   : str
            default total to perform the analysis on total sensitivity index.
            give any other string to do it on first order index
    number  : int
            find the n number of reactions with highest sensitivity indices
    Returns
    -----------
    df      : dataframe
            a dataframe with reaction number and sensitivity indices
    """

    sens_dfs = parse_sa_analysis(path)

    df = sens_dfs[species][0]
    # df = sens_dfs['glycoaldehyde'][0]
    top = int(number)
    # Initialize boolean checks and check dataframe structure
    if (('S1' not in df) or ('ST' not in df) or ('Parameter' not in df) or
            ('ST_conf' not in df) or ('S1_conf' not in df)):
        raise Exception('Dataframe not formatted correctly')

    if coeff == 'total':
        df['placeholder'] = abs(df.ST)
    else:
        df['placeholder'] = abs(df.S1)

    df = df.nlargest(top, 'placeholder')
    df = df.drop('placeholder', axis=1)
    df = df.reset_index(drop=True)
    return df


def build_reaction(reac):
    """
    This function gets the reactants and products in smiles format
    from reaction list and writes the reaction in a way that
    can be fed to rdkit.
    Parameters
    ------
    reac    : str
            A reaction string with reactants/products and their
            stoichiometric coeffs. Typically an entry from complete
            reaction list
    Returns
    --------------
    reaction_list   : str
                    list of reactants and products of a specific reaction
                    seperated by space
    """

    reactants = []
    products = []
    for species in reac:
        if (float(species[0])) < 0:
            reactants.append(species[1])
        else:
            products.append(species[1])

    seperator = ['']
    reaction_list = reactants + seperator + products

    return reaction_list


def draw_top_reactions(path, species, complete_list,
                       number=5, save_image='yes'):
    """
    This function draws the top n reactions with highest sensitivity indices
    with the help of rdkit so that the user can visualize what are the most
    sensitive reactions.
    Parameters
    ------------
    path        : str
                the file path where all the files generated from
                sensitivity analysis are stored.
    species     : str
                name of the species or output measure on which sensitivy
                analysis will be performed
    complete_list    : list
                      A list of all the reactions with reactant and
                      product species and their stoichiometric coeffs
    number      : int
                now many reactions to draw, default = 5
    save_image  : boolean
                to save the image in a .png file
    Returns
    ------------------
    A png file with all the reactions
    """
    df_top = get_top_ones(path, species, number=number)
    parameters = df_top.Parameter.values
    reaction_numbers = [int(parameters[i][1:]) for i in range(number)]

    if save_image:
        save_top_reactions(path, reaction_numbers, species, complete_list)

    # plot_sensitivity_results(df_top)


#        print(reaction)
#    return rxn_strings


def save_top_reactions(image_path, reactions, species, complete_list):
    """
    This function save the reactions drawn using rdkit in a png file
    Parameters
    -------------
    image_path  : str
                the file path where all the files generated from
                sensitivity analysis are stored.
    reactions   : list
                A list of reaction numbers to be drawn
    species     : str
                name of the species or output measure on which sensitivy
                analysis will be performed
    complete_list    : list
                      A list of all the reactions with reactant and
                      product species and their stoichiometric coeffs
    Returns
    ------------------
    A png file with all the reactions
    """
    rxn_strings = []

    for number in reactions:
        reaction = complete_list[number]
        rxn_strings.append(build_reaction(reaction))

    #    rxns = [ChemReac.ReactionFromSmarts(rxn_str) for rxn_str in rxn_strings]

    #    ind_image = ReactionToImage(rxns[1])
    ind_image = [reaction_to_image(rxn) for rxn in rxn_strings]

    rxn_image = stitch_pils_vertically(ind_image)
    rxn_image.save(image_path + 'sensitive_reactions_{}.png'.format(species))


def plot_sensitivity_results(path, species, number=5, coeff='total'):
    """
    Draw a bar plot for the reactiions with highest sensitivity indices.
    Parameters
    --------------
    path        : str
                the file path where all the files generated from
                sensitivity analysis are stored.
    species     : str
                name of the species or output measure on which sensitivy
                analysis will be performed
    number      : int
                now many reactions to draw, default = 5
    coeff       : str
                plot sensitivity indices calculated using total
                sensitivity analysis. Default is total, use any other
                string to do it foe first order
    Returns
    ---------------
    A bar plot in the screen.
    """
    df = get_top_ones(path, species, number=number)
    parameter = list(df.Parameter)
    if coeff == 'total':
        sensitivity = list(df.ST)
#   order = np.array(['ST'] * len(df)
        confidence = list(df.ST_conf)
    else:
        sensitivity = list(df.S1)
        #   order = np.array(['ST'] * len(df)
        confidence = list(df.S1_conf)
    #            yerr = [x + e for x,e in zip(sensitivity, confidence) ]
    #            lower = [x - e/2 for x,e in zip(sensitivity, confidence) ]

    barwidth = 0.8

    plt.bar(parameter, sensitivity, width=barwidth, color="#962980",
            edgecolor='black', yerr=confidence, capsize=7)
    plt.ylabel('Sensitivity Index')
    plt.xlabel('Forward Reaction Rates')

