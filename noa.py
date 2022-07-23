import json
from nearly_orthogonal_array.lib_oa import generate, columnRename
# from utilities import configRead


def configCheck(config):
    """check the config"""
    return True


def get_oa(config):
    """get the OA/NOA design"""
    variables = config["variables"]
    N_runs = config["N_runs"]
    levels = []
    varLists = []
    for var in variables:
        if variables[var]["type"] == "FLOAT":
            levels.append(variables[var]["level"])
            varLists.append(var)
        elif variables[var]["type"] == "INT":
            levels.append(variables[var]["level"])
            varLists.append(var)
        elif variables[var]["type"] == "CATEGORY":
            levels.append(len(variables[var]["list"]))
            varLists.append(var)
        elif variables[var]["type"] == "NOMINAL":
            levels.append(len(variables[var]["list"]))
            varLists.append(var)
    design = generate(levels, N_runs, verbose=0)
    design = columnRename(design)

    return design, varLists


def oa2levels(design, varLists, config):
    """transfer OA/NOA design to levels"""
    NOA_exp = {}
    variables = config["variables"]
    # # print(variables)
    N_runs = config["N_runs"]
    design = design
    varLists = varLists

    for i in range(0, len(varLists), 1):
        if variables[varLists[i]]["type"] == "FLOAT" and config["TouchBoundary"] == "True":
            # print(varLists[i])
            noa_levels = (design.T[i]) * (variables[varLists[i]]["max"] - variables[varLists[i]]["min"]) / (
                variables[varLists[i]]["level"] - 1) + [variables[varLists[i]]["min"]] * N_runs
            NOA_exp[varLists[i]] = list(noa_levels)
        elif variables[varLists[i]]["type"] == "FLOAT" and config["TouchBoundary"] == "False":
            # print(varLists[i])
            noa_levels = (design.T[i] + [1] * N_runs) * (variables[varLists[i]]["max"] - variables[varLists[i]]
                                                         ["min"]) / (variables[varLists[i]]["level"] + 1) + [variables[varLists[i]]["min"]] * N_runs
            NOA_exp[varLists[i]] = list(noa_levels)
        elif variables[varLists[i]]["type"] == "INT" and config["TouchBoundary"] == "True":
            # print(varLists[i])
            noa_levels = (design.T[i]) * (variables[varLists[i]]["max"] - variables[varLists[i]]["min"]) / (
                variables[varLists[i]]["level"] - 1) + [variables[varLists[i]]["min"]] * N_runs
            noa_levels = round_list(noa_levels, 0)
            NOA_exp[varLists[i]] = list(noa_levels)
        elif variables[varLists[i]]["type"] == "INT" and config["TouchBoundary"] == "False":
            # print(varLists[i])
            noa_levels = (design.T[i] + [1] * N_runs) * (variables[varLists[i]]["max"] - variables[varLists[i]]
                                                         ["min"]) / (variables[varLists[i]]["level"] + 1) + [variables[varLists[i]]["min"]] * N_runs
            noa_levels = round_list(noa_levels, 0)
            NOA_exp[varLists[i]] = list(noa_levels)
        elif variables[varLists[i]]["type"] == 'CATEGORY' or variables[varLists[i]]["type"] == 'NOMINAL':
            # print(varLists[i])
            noa_category_levels = []
            for j in range(0, N_runs, 1):
                noa_category_levels.append(
                    variables[varLists[i]]["list"][design.T[i][j]])
            NOA_exp[varLists[i]] = list(noa_category_levels)
    return NOA_exp


def round_list(l, d):
    """round all list"""
    return [int(round(i, d)) for i in l]


def get_noa_exp(config):
    """get nearly orthogonal array exp."""
    design, varLists = get_oa(config)
    # print(design.T)
    # print(varLists)
    noa_exp = oa2levels(design, varLists, config)
    # print(NOA_exp)
    return noa_exp


# config = configRead('config/OA_config.json', configCheck)
# get_noa_exp(config)
