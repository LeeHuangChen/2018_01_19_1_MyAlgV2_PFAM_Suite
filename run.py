import configurations as conf
from src import GenFasta, util
from src import blast_suite as blast
from src import build_HSPInt_graph as buildGraph
from src import defineBordersFromGraph as findBorders
from src import visualization as vis
from src import PFAMComparison as pfamComp
from cPickle import dump
import os
import datetime


def read_families():
    families = []
    with open(conf.familiesFile, "r") as f:
        for i, line in enumerate(f):
            # format:
            # 2840	C1-set
            arr = line.split("\t")
            numprot = int(arr[0].strip())
            fam_name = arr[1].strip()
            if numprot > 50:
                families.append((numprot, fam_name))
    return families


def runAlg(FamNames, filename):
    # fam_name=family[1]
    # numprot=family[0]

    # generate fasta sequence
    outfolder = conf.fastaFolder
    util.generateDirectories(outfolder)
    #filename = str(numprot) + "_" + fam_name

    outdir = os.path.join(outfolder, filename)
    GenFasta.GenerateFastaInputForMultiFamilies(FamNames, outdir)

    # generate protein lengths
    plenFolder = conf.proteinLenFolder
    util.generateDirectories(plenFolder)
    plenDict = blast.generateProtLenDict(conf.fastaFolder, filename)
    dump(plenDict, open(os.path.join(plenFolder, filename), "wb"))

    # create blast databases
    blast.makeblastdb(conf.fastaFolder, filename)

    # conduct all to all BLASTp
    alltoallFolder = conf.alltoallFolder
    util.generateDirectories(alltoallFolder)
    blast.alltoallBlastP(conf.fastaFolder, filename,os.path.join(alltoallFolder, filename))

    # build HSPIntGraph
    seqSimGraph = buildGraph.build_graph(filename, conf.alltoallFolder)

    # identify protein module borders
    # putative domains
    numModules, moduleFamilyInfo = findBorders.generatePutativeModules(seqSimGraph)
    putativeResult = vis.visualizeModuleFamilyInfo(moduleFamilyInfo)

    # remove submodules
    findBorders.removeSubModules(moduleFamilyInfo)
    moduleResult = vis.visualizeModuleFamilyInfo(moduleFamilyInfo)
    # print moduleResult

    # rename modules to have lower numbers
    numModulesAfterFilter = findBorders.renameModules(moduleFamilyInfo)
    moduleResultRenamed = vis.visualizeModuleFamilyInfo(moduleFamilyInfo)
    # print "numberOfModules: ", numModulesAfterFilter
    # print moduleResultRenamed

    # output the results
    util.generateDirectories(conf.resultsFolder)

    consizePath = os.path.join(conf.resultsFolder, filename + "_Modules.txt")
    with open(consizePath, "w") as f:
        f.write(moduleResultRenamed)

    detailedPath = os.path.join(conf.resultsFolder, filename + "_detailedResults.txt")
    with open(detailedPath, "w") as f:
        f.write("Putative Domains: " + str(numModules) + "\n" + putativeResult + "\n")
        f.write("RemoveSubModules: \n" + moduleResult + "\n")
        f.write("Final Module Definition: " + str(numModulesAfterFilter) + "\n" + moduleResultRenamed + "\n")

    # compare the borders with pfam definitions side by side
    comparePFamVis = pfamComp.visualizeResultsWithPFam(moduleFamilyInfo)
    pfamCompPath = os.path.join(conf.resultsFolder, filename + "_pFamSideBySide.txt")
    with open(pfamCompPath, "w") as f:
        f.write(comparePFamVis)


def main():
    # load all families in PFam
    families = read_families()
    families.sort(key=lambda x: x[0])

    # famNames = []
    # util.generateDirectories(conf.resultsFolder)
    #
    # for family in families:
    #
    #     famNames.append(family[1])
    #
    #     filename = datetime.datetime.now().strftime("%Y%m%d_%I%p") + "_" + famNames[0]
    # #runAlg([families[1], families[2]])

    famNames = ["Neur_chan_memb"]
    filename = datetime.datetime.now().strftime("%Y%m%d_%I%p") + "_" + famNames[0]
    runAlg(famNames, filename)


def test():
    pfamComp.generatePFamInfoByProtein()


if __name__ == '__main__':
    # test()
    main()

