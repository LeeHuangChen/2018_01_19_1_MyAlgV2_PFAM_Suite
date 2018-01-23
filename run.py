import configurations as conf
from src import GenFasta, util
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


def main():
    # load all families in PFam
    families = read_families()
    families.sort(key=lambda x: x[0])

    # generate 1 family
    outfolder = conf.fastaFolder
    util.generateDirectories(outfolder)
    fam_name = families[0][1]
    outdir = os.path.join(outfolder, fam_name + ".txt")

    GenFasta.GenerateFastaInputForGivenFamily(fam_name, outdir)

    #


if __name__ == '__main__':
    main()
