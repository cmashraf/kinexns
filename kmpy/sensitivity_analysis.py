import os


def gen_params(sample_number, sa_path, input_file_name, output_file_name):

    input_file = sa_path + input_file_name
    output_file = sa_path + output_file_name
    os.system('python -m SALib.sample.saltelli -n {} -p {} -o {}'.format(sample_number, input_file, output_file))


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
    None
    """
    out_file = path + 'analysis_{}.txt' .format(name)
    params_file = path + 'params.txt'
    result = path + 'model_solutions.txt'
    print(out_file)
    print(processors)
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
