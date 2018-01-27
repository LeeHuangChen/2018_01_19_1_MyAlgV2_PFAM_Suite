import configurations as conf
from src import GenFasta, util
from src import blast_suite as blast
from src import build_HSPInt_graph as buildGraph
from src import defineBordersFromGraph as findBorders
from cPickle import dump
import os


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


def runAlgOnOneFamily(family):
    fam_name=family[1]
    numprot=family[0]

    # generate fasta sequence
    outfolder = conf.fastaFolder
    util.generateDirectories(outfolder)
    famFilename = str(numprot) + "_" + fam_name
    outdir = os.path.join(outfolder, famFilename)
    GenFasta.GenerateFastaInputForGivenFamily(fam_name, outdir)

    # generate protein lengths
    plenFolder = conf.proteinLenFolder
    util.generateDirectories(plenFolder)
    plenDict = blast.generateProtLenDict(conf.fastaFolder, famFilename)
    dump(plenDict, open(os.path.join(plenFolder, famFilename),"wb"))

    # create blast databases
    blast.makeblastdb(conf.fastaFolder, famFilename)

    # conduct all to all BLASTp
    alltoallFolder = conf.alltoallFolder
    util.generateDirectories(alltoallFolder)
    blast.alltoallBlastP(conf.fastaFolder, famFilename,os.path.join(alltoallFolder, famFilename))

    # build HSPIntGraph
    seqSimGraph = buildGraph.build_graph(famFilename, conf.alltoallFolder)

    # identify protein module borders
    findBorders.generatePutativeModules(seqSimGraph)


def main():
    # load all families in PFam
    families = read_families()
    families.sort(key=lambda x: x[0])

    runAlgOnOneFamily(families[0])


if __name__ == '__main__':
    main()
