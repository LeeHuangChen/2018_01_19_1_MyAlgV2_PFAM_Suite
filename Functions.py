import os
import Configurations as conf
import util
from macros import *

def printLine(fileInfo,lineInfo):
	lineIndex=lineInfo[0]
	line=lineInfo[1]
	print line

def printAllLinesInAllFiles():
	forAllLinesInAllFiles(printLine,conf.inputFolder)


def test():
	printAllLinesInAllFiles()
	

def main():
	test()

if __name__ == '__main__':
	main()