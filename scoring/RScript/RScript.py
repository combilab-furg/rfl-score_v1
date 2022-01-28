import os, sys

def get_model_160_path():
    model_path = os.path.join(os.path.join(__file__.rsplit("/", 1)[0], "model"), "scoring_function_160.rda")
    return model_path

def get_model_250_path():
    model_path = os.path.join(os.path.join(__file__.rsplit("/", 1)[0], "model"), "scoring_function_250.rda")
    return model_path

def load_template():
    template_path = os.path.join(os.path.join(__file__.rsplit("/", 1)[0], "templates"), "SF_TEMPLATE.R")
    return template_path

def make_model_script(input_path, output_path, n_features, tmp_folder_name):
    
    sf_template = load_template()
    model_path = None
    
    if n_features == 160:
        model_path = get_model_160_path()
    elif n_features == 250: 
        model_path = get_model_250_path()

    file_path = os.path.join(tmp_folder_name, "SF" + str(n_features) + ".R")
    
    in_file = open(sf_template).read()
    in_file = in_file.replace('#model_path#', model_path)
    in_file = in_file.replace('#tmp_path#', tmp_folder_name)
    in_file = in_file.replace('#input_path#', input_path)
    in_file = in_file.replace('#output_path#', output_path)
    in_file = in_file.replace('#n_features#', str(n_features))

    if os.path.exists(file_path):
        os.remove(file_path)

    out_file = open(file_path, 'w')
    out_file.write(in_file)
    out_file.close()

    return file_path
