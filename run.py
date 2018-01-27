import configurations as conf
from src import GenFasta, util
from src import blast_suite as blast
from src import build_HSPInt_graph as buildGraph
from src import defineBordersFromGraph as findBorders
from src import visualization as vis
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


def runAlg(families):
    # fam_name=family[1]
    # numprot=family[0]

    # generate fasta sequence
    outfolder = conf.fastaFolder
    util.generateDirectories(outfolder)
    #filename = str(numprot) + "_" + fam_name
    filename = datetime.datetime.now().strftime("%Y%m%d_%I%p")+"_"+families[0][1]
    outdir = os.path.join(outfolder, filename)
    GenFasta.GenerateFastaInputForMultiFamilies(families, outdir)

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
    numModules, moduleFamilyInfo = findBorders.generatePutativeModules(seqSimGraph)
    print vis.visualizeModuleFamilyInfo(moduleFamilyInfo)
    findBorders.removeSubModules(moduleFamilyInfo)
    print vis.visualizeModuleFamilyInfo(moduleFamilyInfo)


def main():
    # load all families in PFam
    families = read_families()
    families.sort(key=lambda x: x[0])

    runAlg([families[1], families[2]])


if __name__ == '__main__':
    main()
