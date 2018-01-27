# parameters
# the overlap ratio cutoff to determine if two BLASTp hits on the same protein is considered overlapped
cutoffRatio = .90
evalueCutoff = 1e-4
# The amount of overlap two proteins sequence need to have to be consider the "same" protein sequence
simularProteinRatio = .95

# main core folders
fastaFolder = "Generated/0_FastaSequence"
proteinLenFolder = "Generated/1_ProteinLength"
blastdbFolder = "Generated/2_BlastDB"
alltoallFolder = "Generated/3_BlastAllToAll"
seqSimGraphFolder = "Generated/4_SequenceSimilarityGraphs"

# Sequence Generation
FastaSeqDict = "Resources/FastaSeqDict.cpickle"
FamToArrDictLoc = "Resources/FamToLineDict.cpickle"
familiesFile = "Resources/familyCount.txt"

# Blast all to all
outputformat = 6
alltoallLogFolder = "Generated/log/BlastAllToAll"
blastdbLogFolder = "Generated/log/2_BlastDB"
