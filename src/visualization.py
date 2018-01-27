
def visualizeModuleFamilyInfo(moduleFamilyInfo):
    outStr=""
    for protein in moduleFamilyInfo.keys():
        outStr += protein + ": "
        modules = moduleFamilyInfo[protein]

        for i, module in enumerate(modules):
            m, s, e = str(module[0]), str(module[1]), str(module[2])
            outStr += m + "(" + s + "," + e + ")"
            if i < len(modules)-1:
                outStr += ", "
        outStr += "\n"
    return outStr
