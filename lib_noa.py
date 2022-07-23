"""
Generate orthogonal array (OA) or nearly orthogonal array (NOA).

Method: 1. Check known OA lists
		2. Sequentially adding columns to an existing design with J2-criteria.


Ref: 1. R, SAS, Neilslone
	 2. Xu, H. (2002). An algorithm for constructing orthogonal and nearly-orthogonal arrays with mixed levels and small runs.  Technometrics, 44, 356-368.

"""
import os
import numpy as np
import random
import copy
import math
import csv


def check(levels, weight):
    """return if right inputs"""
    if weight != []:
        if len(levels) != len(weight):
            return False
    return True


def checkCover(Dict1, Dict2):
    """return if Dict1 cover Dict2"""
    flag = True
    for key in Dict2:
        if key not in Dict1:
            flag = False
            break
        elif Dict2[key] > Dict1[key]:
            flag = False
            break
    return flag


def columnRename(design):
    """rename all column"""
    design = design.T
    # print(design.shape)
    for i in range(0, design.shape[0], 1):
        max_i = max(design[i]) + 1
        ind = random.sample(list(np.arange(0, max_i)), max_i)
        # print(ind)
        for j in range(0, design.shape[1], 1):
            design[i][j] = ind[design[i][j]]
    design = design.T
    return design


def deltaMatrix(design, weight):
    """compute delta(i,j)"""
    design = np.matrix(design)
    N_runs = design.shape[1]
    delta_Matrix = np.zeros((N_runs, N_runs))
    for i in range(0, N_runs, 1):
        for j in range(0, N_runs, 1):
            temp = 0
            for k in range(0, design.shape[0], 1):
                if design[k, i] == design[k, j]:
                    temp = temp + weight[k]
            delta_Matrix[i, j] = temp

    return delta_Matrix


def designPrune(designList, target_dict):
    """prune redundant design"""
    new_design = []
    new_levels_exp = []
    design = getDesign(designList["Name"])
    levels_exp = designList["levels_exp"]

    target_dict = copy.copy(target_dict)

    k = 0
    for design_k in design:
        s = levels_exp[k]
        if s in target_dict:
            if target_dict[s] > 0:
                target_dict[s] -= 1
                new_design.append(design_k)
                new_levels_exp.append(s)
        k += 1
    return (new_design, new_levels_exp)


def fillin():
    """fill in all designs"""

    alphabet = ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
                'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    
    s = os.path.dirname(os.path.dirname(__file__))

    csvfile = open(s+'\\OAs\OA_list.csv', 'r', newline='')
    oa_list_file = csv.reader(csvfile)
    header = next(oa_list_file)

    NN = 0

    for line in oa_list_file:
        if os.path.exists(s + '\\OAs\\' + str(line[1]) + '.csv') == 0:
            name = name2levels(line[1])
            nn_level = []
            nn_runs = 1
            for i in range(0, len(name), 1):
                if name[i][1] == 1:
                    nn_level.append(name[i][0])
                    nn_runs = nn_runs * name[i][0]
                else:
                    nn_level = []
            if nn_level != []:
                print(line[1])
                nn_length = 3  # len(name)

                new_file = open(
                    s + '\\OAs\\' + str(line[1]) + '.csv', 'w', newline='')
                w = csv.writer(new_file)
                first_row = []
                for i in range(0, len(name)+1, 1):
                    first_row.append(alphabet[i])
                w.writerow(first_row)
                design = Xu(nn_level, nn_runs, T1=500)
                print(design.T)
                for j in range(0, nn_runs, 1):
                    w.writerow(list([j+1]) +
                               list(design.T[j] + [1] * len(name)))
                    # w.writerow()
                new_file.close
                NN = NN + 1

    csvfile.close
    print(NN)


def findRelevantList(repDict):
    """return lists if cover/under/match"""
    relevantList = {"cover": [], "under": [], "match": None}
    for listDict in OAList:
        thisDict = listDict["levels_dic"]

        coverflag = checkCover(thisDict, repDict)
        underflag = checkCover(repDict, thisDict)

        if coverflag:
            relevantList["cover"].append(listDict)
        if underflag:
            relevantList["under"].append(listDict)
        if coverflag and underflag:
            relevantList["match"] = listDict
    # print(relevantList)
    return relevantList


def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a


def generate(levels, N_runs=None, isOA=False, weight=[], interchangeType="mix", T1=100, T2=10, ratio=1, verbose=1):
    """main process"""
    verbose = verbose
    if not isOA:
        if not check(levels, weight):
            print("Error! The length between levels and weight is different")
            return None
    levels_, levels_dic = levelsRepresented(levels)
    relevantList = findRelevantList(levels_dic)
    levels_exp = None
    design = None

    if isOA:

        # Case 1: N_runs = None, isOA = True
        if not N_runs:
            if verbose:
                print("="*36+"Case 1"+"="*36)
            if relevantList["match"]:
                if verbose:
                    print("There is a match design in database.")
                Name = relevantList["match"]["Name"]
                levels_exp = relevantList["match"]["levels_exp"]
                design = getDesign(Name)
            elif relevantList["cover"]:
                if verbose:
                    print("There are several cover designs in database.")
                coverList = relevantList["cover"]
                min_N_runs = np.inf
                chooseList = None
                for List in coverList:
                    if List["N_runs"] < min_N_runs:
                        chooseList = List
                        min_N_runs = List["N_runs"]

                design, levels_exp = designPrune(chooseList, levels_dic)
            else:
                if verbose:
                    print(
                        "There is not any match/cover design in database, mix process start!")
                result = getMix(levels_dic, verbose=verbose)
                if result is not None:
                    design, levels_exp = result
                else:
                    if verbose:
                        print("Suggest setting isOA = False")
                    return None

        # Case 2: N_runs = #, isOA = True
        else:
            if verbose:
                print("="*36+" Case 2 "+"="*36)
            result = selectMatch(N_runs, levels_dic,
                                 relevantList, verbose=verbose)
            if result is not None:
                design, levels_exp = result
            else:
                if verbose:
                    print(
                        "There is not any design we can offer, try another N_runs or set isOA = False")
                return None

    else:
        if verbose:
            print("="*36+" Case 3 "+"="*36)

        # Case 3: N_runs = None, isOA = False
        if not N_runs:
            N_runs = getLB(levels_dic)

        # Case 4: N_runs = #, isOA = False
        result = selectMatch(N_runs, levels_dic, relevantList, verbose=verbose)
        if result is not None:
            design, levels_exp = result
        else:
            designInit = []
            levelsInit = []
            levelsInit_dic = {}

            max_LB = -np.inf
            if relevantList["under"]:
                if verbose:
                    print("There are several under designs in database.")
                underList = relevantList["under"]
                for List in underList:
                    # print(N_runs, List["N_runs"])
                    if List["N_runs"] == N_runs:
                        Rao_LCM_LB = List["Rao_LCM_LB"]
                        if Rao_LCM_LB > max_LB:
                            max_LB = Rao_LCM_LB
                            designInit = getDesign(List["Name"])
                            levelsInit = List["levels_exp"]
                            levelsInit_dic = List["levels_dic"]
            # print(designInit)
            levels_exp = levelsSeperate(levels, levelsInit, levelsInit_dic)
            design = Xu(levels_exp, N_runs=N_runs, designInit=designInit, weight=weight,
                        interchangeType=interchangeType, T1=T1, T2=T2, ratio=ratio, verbose=verbose)
    design = getOriginDesign(design, levels, levels_exp)
    return np.array(design).T


def getDesign(name):
    """get the design"""

    s = os.path.dirname(os.path.dirname(__file__))
    filename = s + '/OAs/' + str(name) + '.csv'
    csvfile = open(filename, 'r', newline='')
    oa_design_file = csv.reader(csvfile)
    header = next(oa_design_file)

    design = []

    for line in oa_design_file:
        design.append([])
        for i in range(1, len(line), 1):
            design[-1].append(int(line[i])-1)

    csvfile.close
    design = np.array(design)
    """
	# from database

	collect = db['OA_design_R']
	design = []
	cursor = collect.find({"Name":name})
	for line in cursor:
		design = np.array(line["design"])
		# print(design.T)
	"""
    return design.T


def getFullDesign(mixList, dicList):
    """return expanding form in lists"""
    design = []
    levels_exp = []
    count = 0
    for listThis in mixList:
        newDesign = []
        if count == 0:
            newDesign, levels_exp = designPrune(listThis, dicList[count])
            newDesign = copy.deepcopy(newDesign)
            levels_exp = copy.deepcopy(levels_exp)
        else:
            addDesign, addLevels = designPrune(listThis, dicList[count])
            addDesign = copy.deepcopy(addDesign)
            addLevels = copy.deepcopy(addLevels)
            # print(addDesign)
            # print(addLevels)
            m = len(design[0])
            n = len(addDesign[0])
            for i in design:
                newDesign.append(list(i)*n)
            for si in addDesign:
                temp = []
                for l in si:
                    temp += [l]*m
                newDesign.append(temp)
            levels_exp += addLevels
        design = copy.deepcopy(newDesign)
        count += 1
    # print(np.array(design))
    # print(levels_exp)
    return design, levels_exp


def getLB(levels_dic):
    """return rao_lcm_LB"""
    rao_LB = Rao(levels_dic)
    lcm_LB = LCM(levels_dic)
    return int(math.ceil(float(rao_LB)/lcm_LB)*lcm_LB)


def getLeftDic(origin_dic, replace_dic):
    """return left levels_dic after replace"""
    left_dic = {}
    for level in origin_dic:
        x = level
        y = origin_dic[level]
        if level in replace_dic:
            y = max(0, y - replace_dic[level])
        if y > 0:
            left_dic.update({x: y})
    return left_dic


def getLists():
    """get all OA list"""

    s = os.path.dirname(os.path.dirname(__file__))
    csvfile = open(s+'/OAs/OA_list.csv', 'r', newline='')
    oa_list_file = csv.reader(csvfile)
    header = next(oa_list_file)

    OA_list = []

    for line in oa_list_file:
        name = line[1]
        n_runs = int(line[2])
        levels = name2levels(name)
        levels_dic = levels_dic_int(levels)

        OA_list.append({"Name": name, "N_runs": n_runs, "levels": levels, "levels_exp": levels_exp(
            levels), "levels_dic": levels_dic, "Rao_LCM_LB": getLB(levels_dic)})

    csvfile.close
    """
	# from database

	collect = db['OA_list_R']
	OA_list = []
	cursor = collect.find({})

	for line in cursor:
		OA_list.append({"Name": line["Name"], "N_runs": line["N_runs"], "levels": line["levels"], "levels_exp": line["levels_exp"], "levels_dic": line["levels_dic"], "Rao_LCM_LB": line["Rao_LCM_LB"]})
		levels_dic = line["levels_dic"]
		for s in levels_dic:
		    levels_dic[int(s)] = levels_dic.pop(s)
	"""
    return OA_list


def getMax_raoMin_nruns(levels_dic):
    """option 1 of combination"""
    maxRao = -np.inf
    best = None
    temp_raoList = []

    for list_temp in OAList:
        rep_dic = getReplaceDic(levels_dic, list_temp["levels_dic"])
        temp_rao = Rao(rep_dic)
        temp_raoList.append(temp_rao)
        if temp_rao > maxRao:
            maxRao = temp_rao

    minN = np.inf
    for index in range(len(temp_raoList)):
        if temp_raoList[index] == maxRao:
            N = OAList[index]["N_runs"]
            if N < minN:
                minN = N
                best = OAList[index]
    return (getLeftDic(levels_dic, best["levels_dic"]), getReplaceDic(levels_dic, best["levels_dic"]), best, best["N_runs"])


def getMin_heuristicscores(levels_dic):
    """option 2 of combination"""
    minScore = np.inf
    best = None
    for list_temp in OAList:
        score = heuristicScore(
            levels_dic, list_temp["levels_dic"], list_temp["N_runs"])
        if score < minScore:
            minScore = score
            best = list_temp
    return (getLeftDic(levels_dic, best["levels_dic"]), getReplaceDic(levels_dic, best["levels_dic"]), best, best["N_runs"])


def getMix(levels_dic, verbose=1):
    """process of mix combination"""
    temp_dic = copy.copy(levels_dic)
    mixList1 = []
    dicList1 = []
    N1 = 1
    while len(temp_dic):
        temp_dic, rep_dic, listThis, temp_N = getMin_heuristicscores(temp_dic)
        mixList1.append(listThis)
        dicList1.append(rep_dic)
        N1 *= temp_N

    temp_dic = copy.copy(levels_dic)
    mixList2 = []
    dicList2 = []
    N2 = 1
    while len(temp_dic):
        temp_dic, rep_dic, listThis, temp_N = getMax_raoMin_nruns(temp_dic)
        mixList2.append(listThis)
        dicList2.append(rep_dic)
        N2 *= temp_N
    N = None
    mixList = None
    dicList = None
    if N1 < N2:
        mixList = mixList1
        dicList = dicList1
        N = N1
    else:
        mixList = mixList2
        dicList = dicList2
        N = N2
    # print(dicList)
    if verbose:
        print("N in design: {}".format(N))
    nameList = [i["Name"] for i in mixList]
    if verbose:
        print("Combinations in design: {}".format(nameList))
    if N <= 9487:
        return getFullDesign(mixList, dicList)
    else:
        if verbose:
            print(
                "The design is too big(more than 9487) so that we just give you the combinations")
        return None


def getOriginDesign(design, levels, levels_exp):
    """return design with origin order"""
    level_index = []
    designFinal = []
    for level in levels:
        index = -1
        for level_exp in levels_exp:
            index += 1
            if index in level_index:
                continue
            if level_exp == level:
                level_index.append(index)
                break
    for i in level_index:
        designFinal.append(design[i])
    return copy.deepcopy(designFinal)


def getReplaceDic(origin_dic, replace_dic):
    """return replace levels_dic after replace"""
    rep_dic = {}
    for level in origin_dic:
        x = level
        y = origin_dic[level]
        if level in replace_dic:
            y = min(y, replace_dic[level])
            rep_dic.update({x: y})
    return rep_dic


def heuristicScore(levels_dic, levels_dic_list, N_runs):
    """return score of option2"""
    score = N_runs
    for level in levels_dic:
        x = level
        y = levels_dic[level]
        if level in levels_dic_list:
            y = max(0, y - levels_dic_list[level])
        score *= math.pow(x, y)
    return score


def initial(levels, N_runs):
    """generate initial two columns of OA/NOA"""
    ITD = []
    if len(levels) == 1:
        ITD.append([int(i*levels[0]/N_runs) for i in range(N_runs)])
    else:
        ITD_1 = []
        ITD_2 = []
        for i in range(0, N_runs, 1):
            ITD_1.append(int(i*levels[0]/N_runs))
            ITD_2.append(i % levels[1])

        ITD.append(ITD_1)
        ITD.append(ITD_2)

    return np.matrix(ITD)


def interchange(column, design, weight):
    """interchange procedure"""
    J = 0
    time_t = 0
    length = design.shape[1]

    d_M = deltaMatrix(design, weight)

    temp_max = 1
    while time_t < 20 and temp_max > 0:
        d_c = deltaMatrix(column, weight)
        temp_a = -1
        temp_b = -1
        temp_max = 0

        for a in range(0, length-1, 1):
            for b in range(a+1, length, 1):
                if column[a] != column[b]:
                    temp = 0
                    for j in range(0, length, 1):
                        if (j != a) and (j != b):
                            temp = temp + \
                                (d_M[a, j] - d_M[b, j]) * \
                                (d_c[a, j] - d_c[b, j])

                    if temp > temp_max:
                        temp_max = temp
                        temp_a = a
                        temp_b = b
        # print(temp_a, temp_b, temp_max * 2 * weight[design.shape[0]])
        if temp_max > 0:
            temp_c = column[temp_a]
            column[temp_a] = column[temp_b]
            column[temp_b] = temp_c
        # print(J2d(np.vstack((design,column)), weight), time_t)

        time_t = time_t + 1
    J = J2d(np.vstack((design, column)), weight)
    return J, column


"""unnecessary"""


def interchange_bubble(column, design, weight):
    """interchange procedure"""
    J = 0
    time_t = 0
    length = design.shape[1]

    d_M = deltaMatrix(design, weight)

    temp_max = 1
    while time_t < 20 and temp_max > 0:
        d_c = deltaMatrix(column, weight)
        temp_max = 0

        for a in range(0, length-1, 1):
            for b in range(a+1, length, 1):
                if column[a] != column[b]:

                    temp = 0
                    for j in range(0, length, 1):
                        if (j != a) and (j != b):
                            temp = temp + \
                                (d_M[a, j] - d_M[b, j]) * \
                                (d_c[a, j] - d_c[b, j])
                    # print(a, b, temp)
                    if temp > temp_max:
                        temp_max = temp
                        c = column[a]
                        column[a] = column[b]
                        column[b] = c

        time_t = time_t + 1
    J = J2d(np.vstack((design, column)), weight)
    return J, column, time_t


"""unnecessary"""


def interchange_bubble2(column, design, weight):
    """interchange procedure"""
    J = 0
    time_t = 0
    length = design.shape[1]

    d_M = deltaMatrix(design, weight)

    temp_max = 1
    while time_t < 200 and temp_max > 0:
        d_c = deltaMatrix(column, weight)
        temp_max = 0

        skipList = []
        for a in range(0, length-1, 1):
            if a in skipList:
                continue
            for b in range(a+1, length, 1):
                if column[a] != column[b]:

                    temp = 0
                    for j in range(0, length, 1):
                        if (j != a) and (j != b):
                            temp = temp + \
                                (d_M[a, j] - d_M[b, j]) * \
                                (d_c[a, j] - d_c[b, j])
                    # print(a, b, temp)
                    if temp > temp_max:
                        temp_max = temp
                        skipList.append(b)
                        c = column[a]
                        column[a] = column[b]
                        column[b] = c

        time_t = time_t + 1
    J = J2d(np.vstack((design, column)), weight)
    return J, column, time_t


"""unnecessary"""


def interchange_mix(column, design, weight):
    """interchange procedure"""
    J = 0
    time_t = 0
    length = design.shape[1]

    d_M = deltaMatrix(design, weight)

    temp_max = 1
    while time_t < 2 and temp_max > 0:
        d_c = deltaMatrix(column, weight)
        temp_a = -1
        temp_b = -1
        temp_max = 0

        for a in range(0, length-1, 1):
            for b in range(a+1, length, 1):
                if column[a] != column[b]:
                    temp = 0
                    for j in range(0, length, 1):
                        if (j != a) and (j != b):
                            temp = temp + \
                                (d_M[a, j] - d_M[b, j]) * \
                                (d_c[a, j] - d_c[b, j])
                    if temp > temp_max:
                        temp_max = temp
                        temp_a = a
                        temp_b = b
        # print(temp_a, temp_b, temp_max * 2 * weight[design.shape[0]])
        if temp_max > 0:
            temp_c = column[temp_a]
            column[temp_a] = column[temp_b]
            column[temp_b] = temp_c
        # print(J2d(np.vstack((design,column)), weight), time_t)
        time_t = time_t + 1

    while time_t < 200 and temp_max > 0:
        d_c = deltaMatrix(column, weight)
        temp_max = 0

        skipList = []
        for a in range(0, length-1, 1):
            if a in skipList:
                continue
            for b in range(a+1, length, 1):
                if column[a] != column[b]:

                    temp = 0
                    for j in range(0, length, 1):
                        if (j != a) and (j != b):
                            temp = temp + \
                                (d_M[a, j] - d_M[b, j]) * \
                                (d_c[a, j] - d_c[b, j])
                    # print(a, b, temp)
                    if temp > temp_max:
                        temp_max = temp
                        skipList.append(b)
                        c = column[a]
                        column[a] = column[b]
                        column[b] = c

        time_t = time_t + 1
    J = J2d(np.vstack((design, column)), weight)
    return J, column, time_t


def interchange_mix2(column, design, weight):
    """interchange procedure"""
    J = 0
    time_t = 0
    length = design.shape[1]

    d_M = deltaMatrix(design, weight)

    temp_max = 1
    while time_t < 20 and temp_max > 2:
        d_c = deltaMatrix(column, weight)
        temp_a = -1
        temp_b = -1
        temp_max = 0

        for a in range(0, length-1, 1):
            for b in range(a+1, length, 1):
                if column[a] != column[b]:
                    temp = 0
                    for j in range(0, length, 1):
                        if (j != a) and (j != b):
                            temp = temp + \
                                (d_M[a, j] - d_M[b, j]) * \
                                (d_c[a, j] - d_c[b, j])
                    # print(a, b, temp)

                    if temp > temp_max:
                        temp_max = temp
                        temp_a = a
                        temp_b = b
        # print(temp_a, temp_b, temp_max * 2 * weight[design.shape[0]])
        if temp_max > 0:
            temp_c = column[temp_a]
            column[temp_a] = column[temp_b]
            column[temp_b] = temp_c
        # print(J2d(np.vstack((design,column)), weight), time_t)
        time_t = time_t + 1

    while time_t < 200 and temp_max > 0:
        d_c = deltaMatrix(column, weight)
        temp_max = 0

        skipList = []
        for a in range(0, length-1, 1):
            if a in skipList:
                continue
            for b in range(a+1, length, 1):
                if column[a] != column[b]:

                    temp = 0
                    for j in range(0, length, 1):
                        if (j != a) and (j != b):
                            temp = temp + \
                                (d_M[a, j] - d_M[b, j]) * \
                                (d_c[a, j] - d_c[b, j])
                    # print(a, b, temp)
                    if temp > temp_max:
                        temp_max = temp
                        skipList.append(b)
                        c = column[a]
                        column[a] = column[b]
                        column[b] = c

        time_t = time_t + 1
    J = J2d(np.vstack((design, column)), weight)
    return J, column


def J2d(design, weight=[]):
    """compute J2(d) of a design d"""
    if weight == []:
        weight = np.ones(design.shape[0])
    delta_Matrix = deltaMatrix(design, weight)
    J_value = 0
    for i in range(0, delta_Matrix.shape[0], 1):
        for j in range(0, delta_Matrix.shape[1], 1):
            if i < j:
                J_value = J_value + delta_Matrix[i, j] ** 2

    return J_value


def LCM(levels_dic):
    """return LCM lower bound"""
    lcm_temp = 1
    s_bag = []
    for s in levels_dic:
        temp = s
        if levels_dic[s] > 1:
            lcm_temp = lcm(lcm_temp, s*s)
        for sb in s_bag:
            lcm_temp = lcm(lcm_temp, sb*s)
        s_bag.append(s)
    # print(s_bag)
    return lcm_temp


def lcm(a, b):
    """return lowest common multiple."""
    return a * b // gcd(a, b)


def lcmm(*args):
    """Return lcm of args."""
    return reduce(lcm, args)


def levels_dic(levels):
    dic = {}
    for i in range(0, len(levels), 1):
        dic.update({str(levels[i][0]): levels[i][1]})
    return dic


def levels_dic_int(levels):
    dic = {}
    for i in range(0, len(levels), 1):
        dic.update({levels[i][0]: levels[i][1]})
    return dic


def levels_exp(levels):
    exp = []
    for i in range(0, len(levels), 1):
        for j in range(0, levels[i][1], 1):
            exp.append(levels[i][0])
    return exp


def levelsRepresented(levels):
    """return levels with different data structure"""
    levels_dic = {}
    for level in levels:
        if level not in levels_dic:
            levels_dic.update({level: 1})
        else:
            levels_dic[level] += 1
    # print(levels_dic)
    levels_ = []
    for key in sorted(levels_dic, reverse=True):
        levels_.append([key, levels_dic[key]])
    # print(levels_dic)
    return levels_, levels_dic


def levelsSeperate(levels, levels_under, levels_dict):
    """return remaining columns to be design"""
    levels_sep = copy.copy(levels_under)
    levels_dict = copy.copy(levels_dict)
    levels = sorted(levels, reverse=True)
    for level in levels:
        if level in levels_dict:
            if levels_dict[level] > 0:
                levels_dict[level] -= 1
                continue
        levels_sep.append(level)
    return levels_sep


def J2d_LB(levels, N_runs, weight=[]):
    """compute the lower bound L(k) of J2(d) when k = 1, 2, ..., n"""
    levels = np.array(levels)
    n_factors = levels.size
    if weight == []:
        weight = np.ones(n_factors)
    L = np.zeros(levels.size)

    L_1 = 0
    L_2 = 0
    L_3 = 0

    for i in range(0, levels.size, 1):

        L_1 = L_1 + float(N_runs)/levels[i] * weight[i]
        L_2 = L_2 + (levels[i]-1) * (float(N_runs)/levels[i] * weight[i]) ** 2
        L_3 = L_3 + weight[i]

        L[i] = (L_1 ** 2 + L_2 - N_runs * L_3 ** 2)/2

    return L


def name2levels(Name):
    levels = []
    name_list = str(Name).split('.')

    for i in range(1, len(name_list), 2):
        levels.append([])
        levels[-1].append(int(name_list[i]))
        levels[-1].append(int(name_list[i+1]))

    return levels


def newColumn(level_k, N_runs):
    """generate a new balanced column"""
    nC = []
    for i in range(0, N_runs, 1):
        nC.append(int(i*level_k/N_runs))
        # nC.append(i%level_k)
    random.shuffle(nC)
    return nC


"""unnecessary"""


def randomGenerate(levels, N_runs, T=10, ratio=1, weight=[]):
    """generate OA/NOA randomly"""
    levels = np.array(levels)
    n_factors = levels.size
    if weight == []:
        weight = np.ones(n_factors)

    LB = J2d_LB(levels, N_runs, weight)
    # print(LB)
    design = initial(levels, N_runs)

    print("")
    print("======== Initial two columns ========")
    print("")
    print(LB[0], "<<== Lower Bound")
    print(J2d(design[0], weight), design[0])
    print(LB[1], "<<== Lower Bound")
    print(J2d(design, weight), design[1])
    print("")
    print("======== Adding new columns ========")
    print("")

    for k in range(2, n_factors, 1):

        temp_J = np.inf
        print(LB[k], "<<== Lower Bound")
        temp_column = []

        T_T = T

        while T_T >= 1 and temp_J/LB[k] > ratio:
            new_C = newColumn(levels[k], N_runs)
            J = J2d(np.vstack((design, new_C)), weight)
            if J < temp_J:
                temp_J = J
                temp_column = new_C
            T_T = T_T - 1
        print(temp_J, temp_column)
        design = np.vstack((design, temp_column))
        # print(design)

    print("")
    print("======== Done ========")

    J2d_value = J2d(design, weight)
    design = np.matrix(design)
    return J2d_value, design.T
    # return design


def Rao(levels_dic):
    """return RAO lower bound"""
    rao = 1
    for level_dic in levels_dic:
        rao += (level_dic-1)*levels_dic[level_dic]
    return rao


def selectMatch(N_runs, levels_dic, relevantList, verbose=1):
    """return match/cover in database"""
    if relevantList["match"]:
        if relevantList["match"]["N_runs"] == N_runs:
            if verbose:
                print("There is a match design in database.")
            Name = relevantList["match"]["Name"]
            return (getDesign(Name), relevantList["match"]["levels_exp"])
    if relevantList["cover"]:
        if verbose:
            print("There are several cover designs in database.")
        coverList = relevantList["cover"]
        for List in coverList:
            if List["N_runs"] == N_runs:
                return designPrune(List, levels_dic)
    return None


def Xu(levels, N_runs, designInit=[], weight=[], interchangeType="mix", T1=100, T2=10, ratio=1, verbose=1):
    """generate OA/NOA"""
    levels = np.array(levels)
    n_factors = levels.size
    if weight == []:
        weight = np.ones(n_factors)
    if n_factors > 87 and verbose:
        print("="*72)
        print("Warning! The length of factors(levels) is so big(more than 87) that take long time to compute")
    if N_runs > 50 and verbose:
        print("="*72)
        print("Warning! The length of N_runs is so big(more than 50) that take long time to compute")

    LB = J2d_LB(levels, N_runs, weight)
    # print(LB)

    if designInit == []:
        design = initial(levels, N_runs)
    else:
        if verbose:
            print("there is a good init!")
        design = np.matrix(designInit)
    init_k = len(design)
    # print(design)

    if verbose:
        print("")
        print("======== Initial ========")
        print("")
        for i in range(init_k):
            print(LB[i], "<<== Lower Bound")
            print(J2d(design[:i+1], weight), design[i])
        print("")
        print("======== Adding new columns ========")
        print("")

    flag = True
    for k in range(init_k, n_factors, 1):

        temp_J = np.inf
        if verbose:
            print(LB[k], "<<== Lower Bound")
        temp_column = []

        if flag:
            T_T = T1
        else:
            T_T = T2

        new_C = newColumn(levels[k], N_runs)
        while T_T >= 1 and abs(temp_J/LB[k]) > ratio:
            random.shuffle(new_C)

            if interchangeType == "origin":
                (J, col) = interchange(new_C, design, weight)
            else:
                (J, col) = interchange_mix2(new_C, design, weight)

            if J < temp_J:
                temp_J = J
                temp_column = col

            T_T = T_T - 1
        if verbose:
            print(temp_J, temp_column)
        design = np.vstack((design, temp_column))
        # print(design)
        if temp_J > LB[k]:
            flag = False
    if verbose:
        print("")
        print("======== Done ========")
    # print(record)

    J2d_value = J2d(design, weight)
    design = np.array(design)
    # return J2d_value, design.T
    return design


OAList = getLists()

# ================================================================================
# levels_, levels_exp, levels_dic = levelsRepresented([2,2,2,2,4])

# designFinal = getOriginDesign([[1],[2],[3],[4],[5]],[2,2,4,2,2],[2,2,2,2,4])
# designFinal = getOriginDesign([[1],[2],[3],[4],[5],[6]],[2,2,3,3,4,4],[4,3,4,2,3,2])
# print(designFinal)

# levels_sep = levelsSeperate([5,32,1,5,63,2],[5,2,1],{5:1,2:1,1:1})
# design = LCM({4:1,2:1,3:1})

# design = generate([2]*20, isOA = True)
# design = generate([2]*20+[5]*6+[3]*15, isOA = True)
# design = generate([2]*2+[5]*3+[3]*3, isOA = True)
# design = generate([2]*3+[3]*2+[5]*6+[8], isOA = True)
# design = generate([2,3,2,3,2,2,2,2,2,2,2,3,3,2,3,5,7,5,5]+[3]*50+[7]*8, isOA = True)
# design = generate([2]*20+[3]*15+[4]*10+[5]*8+[6]*4+[7]*2, isOA = True)
# design = generate([2,2,4,2,2], N_runs=8, isOA = True)
# design = generate([6]*2+[3]+[5]*9, T1=5, T2=2)
# design = generate([2]*6+[3])
# design = generate([6,7])
# design = Xu([6,7], 42)
# design = generate([3]*14, N_runs=27)
# design = generate([2]*5, N_runs=51)

# design = Xu([2,4], N_runs=16)
# design = Xu([2,2,2,2,4,5], N_runs=8, designInit=getDesign("L8.2.4.4.1"))

# design = designPrune(OAList[1],{4:1,2:1})
# design = getLB({4:1,2:1})
# print(design)
# print(columnRename(design))
# print(getLB({6:2,3:1,5:9}))

# print(getLeftDic({2:6,3:10,5:1},{2:4,5:2}))
# print(getReplaceDic({2:6,3:10,5:1},{2:4,5:2}))

# getBest({2:4,3:7,5:9})
# getMaxMin({2:20,3:15,4:10,5:8,6:4,7:2})
# getMix({2:20,3:15,4:10,5:8,6:4,7:2})
# design, levels_exp = getMix({2:20,3:15,5:6})
# print(np.array(design).shape)
# print(levels_exp)
# fillin()
# getFullDesign([OAList[0],OAList[0]], [{2:2},{2:2}])

# (J,a) = generate( [2,2,2,2,3], 12)
# (J,a) = generateOA.generate( [3,3,3,3,3,3,3,2], 18, 50)
# (J,a) = generateOA.generate( [5,5,5,5,5,5,5,5,5], 50, 100)
# (J,a) = generateOA.generate( [5,5,5,5,5,5], 25, T1=100)
# (J,a) = generateOA.generate( [3,5,2,3], 15 )
# (J,a) = generateOA.generate( [5,2,2,2,2,2,2,2,2], 20, 100)
# # b = generateOA.J2d_LB( [2,2,3,3,3], 18, [1,1,1,1,1])
# b = generateOA.J2d_LB([3,2,2,2,2,2,2,2,2,2], 12, [1,1,1,1,1,1,1,1,1,1])
# # b = generateOA.J2d_LB([3,2,2,2,2], 12, [1,1,1,1,1])
# b = generateOA.J2d_LB([3,5,2,3], 6, [1,1,1,1])
# c = generateOA.initial([3,2,2], 12)
# c = initial([3,5,2,3], 6)
# d = newColumn(2, 12)
# e = deltaMatrix(c, [1,1,1,1])
# f = J2d(c)
# (g1,g2) = interchange(d, c, [1,1,1,1,1])
# h = deltaMatrix(d, [1])
# m = randomgenerateOA.generate([3,3,3,3,3,3,3,2], 18, 500)

# print(a)
# print(b)
# print(c)
# print(d)
# print(e)
# print(f)
# print(g1,g2)
# print(h)
# print(m)


# # exp
# lev = [[3]*7+[2],[6]+[3]*6,[5]*6,[3]*13,[2]*23]
# NN_runs = [18,18,25,27,24]
# TT_time = 100

# JJ = 0.0
# S_N = 0

# test_index = 4
# test_number = 100

# L_B = J2d_LB(lev[test_index], NN_runs[test_index])[-1]

# for nn in range(0, test_number, 1):
# 	(J,a) = generate(lev[test_index], NN_runs[test_index], T1 = TT_time)
# 	# (J,a) = randomgenerate(lev, NN_runs, TT_time)
# 	JJ = JJ + J
# 	if J == L_B:
# 		S_N = S_N + 1

# # print(L_B)
# print(float(S_N)/test_number)
# print(JJ/test_number/L_B)

# # print('============================================')

# End = time.time()
# Elapsed = End - Start

# print(('Operation Time : ' + str(Elapsed)))
# # print("test")

# print(getLists())
# print(getDesign("L32.2.3.4.7.8.1"))
