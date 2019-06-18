import os


def gen_params(sample_number, sa_path, input_file_name, output_file_name):

    input_file = sa_path + input_file_name
    output_file = sa_path + output_file_name
    os.system('python -m SALib.sample.saltelli -n {} -p {} -o {}'.format(sample_number, input_file, output_file))


