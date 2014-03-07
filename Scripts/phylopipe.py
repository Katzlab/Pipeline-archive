#!/usr/bin/python
from Pipeline import Pipeline

PathtoFiles ='../Files/'


def main(): 	
	print '################################################################################'
	print 'This script assumes your data files and scripts are in the folders they came in. '
	print 'It also assumes there is a list of OGs you are interested in in the Files folder.' 
	print '################################################################################'
	testPipelineList = raw_input('What is the name of your OGofInterestList? ')
	'''
	#Can be reactivated to allow restarting
	restart  = raw_input('Are you restarting an analysis? (y/n) ')
	if restart[0] == 'y':
		parallel = raw_input('Would you like to run the file (manually) in parallel? (y/n) ')
		if parallel[0] == 'y':
			newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, 'parallel')	
		else:
			newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, 'yes')		
	else:
		newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, 'no')	
	'''
	newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, 'no')	
main()
