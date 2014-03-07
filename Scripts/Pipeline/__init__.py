import os
from Taxon import Taxon
from Gene import Gene


##############################################################

class Pipeline: #given a file with a list of OGs, makes the interestList and
					#a fasta file for the OGs.  Also a series of paths for output
					#and a list of Taxon objects
	
	def __init__(self,OGofInterestFile, PathtoFiles, restart):
		self.PathtoFiles = PathtoFiles
		self.OGofInterestFile = OGofInterestFile
		self.OGList = [] #list of gene names
		self.OGsofInterestList = [] #list of Gene instances
		self.TaxonList = [] # list of Taxon instances
		self.OGswithTaxofInterest = []
		#make directories and set paths:
		PATH = '../' + OGofInterestFile.split('/')[-1] + '/'
		os.system('mkdir ' +  PATH)
		os.system('mkdir ' + PATH  + '/Temp')
		os.system('mkdir ' + PATH  + '/Output')		
		self.PathtoTemp = PATH + 'Temp'
		self.PathtoOutput = PATH + 'Output'
		self.PathtoSeqFiles = PathtoFiles + '/TaxonDataFiles/'
		self.PathtoBlastFiles = PathtoFiles + '/BlastFiles/'
		self.PathtoOGFiles = PathtoFiles + '/allOG5Files/'
		
		
		if restart == 'no':
			self.setInterestList()
			self.writelog('Pipeline:' + OGofInterestFile + ',' + PathtoFiles)
			self.makedb()
			self.populate()
		else:
			if restart == 'add':
				self.setInterestList()
				print 'adding taxa'
				self.addtax(self.PathtoOutput,OGofInterestFile,PathtoFiles)
			elif restart == 'parallel':
				OGswithTaxofInterestin = raw_input('Enter your partial list file name ')
				self.restartP(self.PathtoOutput,OGofInterestFile,PathtoFiles,OGswithTaxofInterestin)
			elif restart == 'test':
				self.test(self.PathtoOutput,OGofInterestFile,PathtoFiles)
			else:
				self.setInterestList()
				self.restart(self.PathtoOutput,OGofInterestFile,PathtoFiles)
		
##############################################################
#runs the pipeline
##############################################################
	def setInterestList(self):	
	
		#to allow running of different Pipelines, output files names have to be different
		#get the list of OGs for this pipeline and build gene instances
		
		for line in open(self.OGofInterestFile,'r'):
			self.OGList.append(line.strip())
			NewGene = Gene(line,self)
			self.OGsofInterestList.append(NewGene)
			
		
	def makedb(self):		
		self.build_orthoMCL_db() #Gets individual fasta files for each OG of interest and also concatenates them into Temp/OrthoMCL_listDB.fas

	def populate(self):	
		i = 0
		#read the taxon files and build Taxon Instances
		for file in os.listdir(self.PathtoSeqFiles):
			#print i
			if file[0] != '.':
				NewTaxon = Taxon(file, self, 'no')
				self.TaxonList.append(NewTaxon)
				self.writelog('Taxon:' + NewTaxon.code)
				for item in NewTaxon.OGlist:
					if item not in self.OGswithTaxofInterest:
						self.OGswithTaxofInterest.append(item)
				#NewTaxon.clearmem()
			i = i + 1	
			
		self.writelog('OGs with Tax: ' + str(self.OGswithTaxofInterest))	
		for gene in self.OGsofInterestList:
			#print gene.OG
			#os.system('rm ' + self.PathtoOutput + '/igp_to_remove.txt_' + gene.OG)
			if gene.OG in self.OGswithTaxofInterest:
				for tax in self.TaxonList:
					gene.setTaxa(tax)
					
			
				gene.getAllSeqs('no')
			
				gene.getSeqCodes()
				gene.runmafft(0) 
				a = gene.run_guidance()
				if a == False:
					self.writelog('Guidance Problem:\nGene:' + gene.OG)
					self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
				else:				
					gene.move_guid()
					b = gene.mask()
					if b == False:
						self.writelog('Guidance Problem:\nGene:' + gene.OG)
						self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
					else:					
						gene.callraxml()
						gene.renameTree()
				
						self.writelog('Gene:' + gene.OG)
				gene.get_sister_tax()
				gene.contamination()
				gene.removeParalogs()
		self.cleanup()
		#pick paralogs here, after checking for contaminants
##############################################################
# methods
##############################################################

	def writelog(self,string):
		self.logfile = open(self.PathtoOutput + '/logfile', 'a') #this will hold info as the pipeline progresses, hopefully to help restart
		self.logfile.write(string + '\n')
		
		
	def build_orthoMCL_db(self):  #Gets individual fasta files for each OG of interest and also concatenates them into Temp/OrthoMCL_listDB.fas
		for OG in self.OGList:
			os.system('cp ' + self.PathtoOGFiles + '/' + OG + ' ' +  self.PathtoTemp )
		os.system('cat ' + self.PathtoTemp + '/OG5_* >> ' + self.PathtoTemp + '/OrthoMCL_listDB.fas')
		

	def cleanup(self):
		os.system('zip -r ' + self.PathtoOutput + '/outputFiles.zip ' + self.PathtoOutput + '/*')
		for folder in ['AllSeqFiles','BestRaxTrees','ContaminationFiles','fasta2keep','ForRaxML','Guidance','logfile','MAFFT','Needlescores','*.fas','*_seqcodes.txt','*_sequencesKept.txt']:
			os.system('rm -r ' + self.PathtoOutput + '/' + folder)
		os.system('rm -r ' + self.PathtoOutput + '/../Temp')
		os.system('rm formatdb.log')
		os.system('rm NoSeqsRemoved')

##############################################################
	#restart
##############################################################
	
	def restart(self,PathtoOutput,OGofInterestFile,pathtoFiles):
		infile = open(PathtoOutput + '/logfile', 'r')
		loglist = infile.readlines()
		lastline = loglist[-1]
		
	
		if lastline.split(':')[0] == 'Pipeline':
			newPipe = Pipeline(OGofInterestFile, pathtoFiles, 'no')	
		if lastline.split(':')[0] == 'Taxon':
			print 'this isnt quite ready -- be prepared to reset once Taxa are done'
			self.restartTaxon(loglist,lastline.split(':')[1])	
		if lastline.split(':')[0] == 'Gene' or lastline.split(':')[0] == 'OGs with Tax':
			self.restartGene(loglist)
			
			
			
	
	def restartTaxon(self,loglist,lastCompleteTaxon):
		taxlist = []
		for line in loglist:
			type = line.split(':')[0]
			if type == 'Taxon':
				taxlist.append(line.split(':')[1].strip())
			
			#NOT FINISHED
			#for tax in taxlist:
			#	newtax =  Taxon(fileName, self, 'yes')
	
	
	
	
	
	
	def restartGene(self,loglist):# GOOD
		genelist = []
		for line in loglist:
			type = line.split(':')[0]
			if type == 'OGs with Tax':
				self.OGswithTaxofInterest = line.split(':')[1]
				#print self.OGsofInterestList
			if type == 'Gene':
				genelist.append(line.split(':')[1].strip())

	
		for gene in self.OGsofInterestList:
			#print gene.OG
			if gene.OG not in genelist:
				if gene.OG in self.OGswithTaxofInterest:
					gene.getAllSeqs('yes')			
					gene.getSeqCodes()
					gene.runmafft(0) #*Pipeline1* Comment out to skip ingroup paralog removal
					a = gene.run_guidance()
					if a == False:
						self.writelog('Guidance Problem:\nGene:' + gene.OG)
						self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
					else:				
						gene.move_guid()
						b = gene.mask()
						if b == False:
							self.writelog('Guidance Problem:\nGene:' + gene.OG)
							self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
						else:					
							gene.callraxml()
							gene.renameTree()
					
							self.writelog('Gene:' + gene.OG)
							
	def restartP(self,PathtoOutput,OGofInterestFile,PathtoFiles,OGswithTaxofInterestin):	
	
		OGswithTaxofInterest = []		
			
		OGswithTaxofInterestfile = open(self.PathtoFiles + '/' + OGswithTaxofInterestin,'r')
		for line in OGswithTaxofInterestfile:
			OGswithTaxofInterest.append(line.strip())
		
		for line in OGswithTaxofInterest:
			print 'line = ' + line
			self.OGList.append(line.strip())
			NewGene = Gene(line,self)
			self.OGsofInterestList.append(NewGene)	
						
		for gene in self.OGsofInterestList:
			#print gene.OG
			#if gene.OG not in genelist:
			print OGswithTaxofInterest
			print gene.OG
			if str(gene.OG) in OGswithTaxofInterest:
				print gene.OG
				gene.getAllSeqs('yes')	
				print 'allSeqs for ' +  gene.OG
				gene.getSeqCodes()
				print 'Seqcodes for ' +  gene.OG
				gene.runmafft(0)
				a = gene.run_guidance()
				if a == False:
					self.writelog('Guidance Problem:\nGene:' + gene.OG)
					self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
				else:				
					gene.move_guid()
					b = gene.mask()
					if b == False:
						self.writelog('Guidance Problem:\nGene:' + gene.OG)
						#self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
					else:					
						gene.callraxml()
						gene.renameTree()
				
						self.writelog('Gene:' + gene.OG)
		
		
		
		
						
##############################################################
#add taxa
##############################################################
	def addtax(self,PathtoOutput,OGofInterestFile,PathtoFiles):
		self.makedb()
		i = 0
		#read the taxon files and build Taxon Instances
		for file in os.listdir(self.PathtoFiles + '/TaxonDataFiles-toAdd'):
			#print i
			if file[0] != '.':
				NewTaxon = Taxon(file, self, 'no')
				self.TaxonList.append(NewTaxon)
				self.writelog('Taxon:' + NewTaxon.code)
				for item in NewTaxon.OGlist:
					if item not in self.OGswithTaxofInterest:
						self.OGswithTaxofInterest.append(item)
			
			i = i + 1	
			
		self.writelog('OGs with Tax Added: ' + str(self.OGswithTaxofInterest))	
		outfile = open(self.PathtoFiles + 'OGs2Redo','a')
		for og in self.OGswithTaxofInterest:
			outfile.write(og + '\n')
			
		outfile.close()
			
		print'make sure files are in place and restart the pipeline in parallel using the file "OGs2Redo"'
		
##############################################################		
	def test_restartwalignment(self,PathtoOutput,OGofInterestFile,PathtoFiles):
		self.setInterestList()
		for gene in self.OGsofInterestList:
				gene.getSeqCodes()
				gene.runmafft(0)
				a = gene.run_guidance()
				if a == False:
					self.writelog('Guidance Problem:\nGene:' + gene.OG)
					self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
				else:				
					gene.move_guid()
					b = gene.mask()
					if b == False:
						self.writelog('Guidance Problem:\nGene:' + gene.OG)
						self.restart(self.PathtoOutput,self.OGofInterestFile,self.PathtoFiles)
					else:					
						gene.callraxml()
						gene.renameTree()
				
						self.writelog('Gene:' + gene.OG)
	def test(self,PathtoOutput,OGofInterestFile,PathtoFiles):
		self.setInterestList()
		for gene in self.OGsofInterestList:
			gene.getseqsfromCodeFile()
			gene.removeParalogs()
			gene.clearmem
##############################################################