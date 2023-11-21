def write_params(
        param_dict,
        param_file
):
    param_f = open(param_file, "w")
    for key in param_dict.keys():
        print("{:<15}{}\n".format(key, param_dict[key]))
        param_f.write("{:<15}{}\n".format(key, param_dict[key]))
    param_f.write("\n\n")
    param_f.close()