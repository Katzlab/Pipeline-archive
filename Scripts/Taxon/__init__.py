import os,re,gc
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Data import CodonTable
from Bio.Data.CodonTable import CodonTable
from Bio import AlignIO 
from Bio import pairwise2
from Bio.Align.Applications import ClustalwCommandline
from Bio.Emboss.Applications import NeedleCommandline
##############################################################

toRemove = []
toKeep = []
blastCutOff = 1e-25
##############################################################

class Taxon: 
	def __init__(self, fileName, Pipeline, restart):
		self.PathtoSeqFiles = Pipeline.PathtoSeqFiles
		self.PathtoBlastFiles = Pipeline.PathtoBlastFiles
		self.PathtoFiles = Pipeline.PathtoFiles
		self.PathtoOutput = Pipeline.PathtoOutput
		self.PathtoTemp = Pipeline.PathtoTemp
		self.OGList  = Pipeline.OGList 
		
		self.fileName = fileName
		self.seqs = [] #a list of seqeunces from the taxon
		self.nuc_prot = "None" #nuc if nucleotide, prot if protein?
		self.blast = [] #a list of blast records		
		self.code = "" #taxon code (get from file?)
		self.OGlist = [] # #list of Ogs of interest in the taxon--'of interest' depends on the pipeline
		self.OGxSeqHash = {} #dictionary with seqeunces for each og of interest are stored
		self.OGxSeqHashProt = {} #dictionary with translated seqeunces for each og of interest are stored
		self.OGxSeqHashtoKeep = {} #dictionary with seqeunces for each og of interest, after ridding, are stored		
		
		
		
		if restart == 'no':
			self.set_seqs(self.PathtoSeqFiles,fileName)
			self.set_code()
			#print self.code
			self.set_blast(fileName, self.PathtoBlastFiles)
			self.parseblast()
		
			if any(a != [] for a in self.OGxSeqHash.values()) == True or any(a != [] for a in self.OGxSeqHashProt.values()) == True:    #  or self.OGxSeqHashProt != []:
				if self.is_nuc() == 'nuc': #need to translate nuc to prot		
					#print 'translating...' + self.code
					self.getbestblastallseqs(1) #translation table = 1.  Someday fix to allow different tables for ciliates, etc.	
				else: 
					for filename in os.listdir(self.PathtoTemp):
						if re.search(self.code, filename):
							os.system('mv ' + self.PathtoTemp + '/' + filename + ' ' + self.PathtoOutput) #move protein files		
				self.pairwiseclustalList() 
			for i in range(len(self.blast)): #tried just setting self.blast = [] but some lists toolong? this removes a big list freeing memory
				self.blast.pop()
			gc.collect()

					
		elif restart == 'YT': #add a method to compare needle scores prir to translating, not after, to keep silent paralogs.
			self.set_seqs(self.PathtoSeqFiles,fileName)
			self.set_code()
			print self.code
			self.set_blast(fileName, self.PathtoBlastFiles)
			self.parseblast()
		
			if any(a != [] for a in self.OGxSeqHash.values()) == True or any(a != [] for a in self.OGxSeqHashProt.values()) == True:    #  or self.OGxSeqHashProt != []:				
				if self.is_nuc() == 'nuc': #need to translate nuc to prot
					self.pairwiseclustalListnt() 
					print 'translating...' + self.code
					self.getbestblastallseqsnt(1) #translation table = 1.  Someday fix to allow different tables for ciliates, etc.	
				else: 
					
					for file in os.listdir(self.PathtoTemp):
						if re.search(self.code, file):
							os.system('mv ' + self.PathtoTemp + '/' + file + ' ' + self.PathtoOutput) #move protein files		
					self.pairwiseclustalList() 
			for i in range(len(self.blast)): #tried just setting self.blast = [] but some lists toolong? this removes a big list freeing memory
				self.blast.pop()
			gc.collect()
			
		elif restart == 'test':
			self.pairwiseclustalListTest()
		
		else:
			infas = SeqIO.parse(open(self.PathtoOutput + '/fasta2keep/' + filename + 'fastatokeep.fas','r'),'fasta')
			if infas != []:
				self.code = infas[0].id.split('_')[0] + '_' + infas[0].id.split('_')[1]
				for item in infas:     
					self.OGxSeqHashtoKeep[(OG,self.code)].append(item)
        
##############################################################
#  methods
##############################################################
		
	def clearmem(self):
		#try to free up memory...			
		self.seqs = [] #a list of seqeunces from the taxon
		self.nuc_prot = "None" #nuc if nucleotide, prot if protein?
		self.blast = [] #a list of blast records		
		self.code = "" #taxon code (get from file?
		self.OGlist = [] # #list of Ogs of interest in the taxon--'of interest' depends on the pipeline
		self.OGxSeqHash = {} #dictionary with seqeunces for each og of interest are stored
		self.OGxSeqHashProt = {} #dictionary with translated seqeunces for each og of interest are stored
		self.OGxSeqHashtoKeep = {} #dictionary with seqeunces for each og of interest, after ridding, are stored		
	
		
		
	def set_seqs(self,PathtoSeqFiles, filename):
		infile = SeqIO.parse(open(PathtoSeqFiles + filename,'r'),'fasta')
		for seq in infile:
			self.seqs.append(seq)	
	
	def set_code(self):
		self.code = self.seqs[0].id.split('_')[0] + '_' + self.seqs[0].id.split('_')[1]+ '_' + self.seqs[0].id.split('_')[2]
		print self.code

		
		
	def is_nuc(self):
		try:
			testSeq = self.seqs[0]

		except:
			self.nuc_prot = "None"
	
		A = testSeq.seq.count('A')
		T = testSeq.seq.count('T')
		G = testSeq.seq.count('G')
		C = testSeq.seq.count('C')

		testSeqCount = float(A+T+C+G)
		length = float(len(testSeq))

		if testSeqCount/length > 0.9:
			self.nuc_prot = 'Nuc'
			return 'nuc'
		else:
			self.nuc_prot = 'Prot'
		

			
	def set_blast(self,fileName, PathtoBlastFiles):	
		test = ''
		for blastFileName in os.listdir(PathtoBlastFiles):
			if re.search(fileName[0:-8], blastFileName):
				test = 'ok'
				#print self.code + ' ok'
				blastrecords = NCBIXML.parse(open(PathtoBlastFiles + blastFileName,'r'))		
				for blast in blastrecords:
						self.blast.append(blast)
		if test == '':
			ans = raw_input('no blast file found for ' + fileName + '. run blast? y/n ')
			if ans == 'y':
				self.get_blast(fileName, PathtoBlastFiles)
				self.set_blast(fileName, PathtoBlastFiles)
			else: 
				self.set_blast(fileName, PathtoBlastFiles)


	def get_blast(self,fileName, PathtoBlastFiles):
		os.system('formatdb -i ' + self.PathtoTemp + '/OrthoMCL_listDB.fas  -p T -o F')
		self.is_nuc()
		#print self.nuc_prot
		if self.nuc_prot == 'Nuc':
			cline = 'blastx -db ' + self.PathtoTemp + '/OrthoMCL_listDB.fas -query ' + self.PathtoSeqFiles+ '/' + fileName + '  -num_descriptions 1 -num_alignments 1 -evalue 1e-15 -outfmt 5 -num_threads 2 -out ' + self.PathtoBlastFiles + '/' + fileName + '_BlastOut'
			os.system(cline)

		elif self.nuc_prot == 'Prot':
			cline = 'blastp -db ' + self.PathtoTemp + '/OrthoMCL_listDB.fas -query ' + self.PathtoSeqFiles+ '/' + fileName + '  -num_descriptions 1 -num_alignments 1 -evalue 1e-15 -outfmt 5 -num_threads 2 -out ' + self.PathtoBlastFiles + '/' + fileName + '_BlastOut'
			os.system(cline)

		
		
		
	def checkFileSize(self):
		count = 0
		infile = open('../Temp/' + arg,'r')
		for line in infile:
			if line[0] == '>':
				count = count + 1
		#print 'count = ' + str(count)
		if count > 50002:
			import s4
			s4.splitf(arg,50000)	
			
			dir = os.pardir + '/Temp'
			for f in os.listdir(dir):
	 			if re.match('\d',f):
					cline = 'blastp -db ../Temp/OrthoMCL_listDB.fas -query ../Temp/' + arg + '  -num_descriptions 1 -num_alignments 1 -evalue 1e-15 -outfmt 5 -num_threads 2 -out ../Temp/Temp/' + arg + '_BlastOut'
					os.system(cline)
					os.system('cat ../Temp/Temp/*_BlastOut > ../Temp/' + arg + '_BlastOut')
					os.system('rm -r ../Temp/Temp')				
		else:
			return True	
			
			
		
###################################################################################################	
#	*Taxon1*: parseblast, getgrpp and makeFasta								  							  #
#	Looks through the blasts and makes files for the taxon for each OG							  #
###################################################################################################	

	def parseblast(self): 

		interestList = self.OGList 
		
		for OGGroup in interestList:
			self.OGxSeqHash[(OGGroup,self.code)] = []
			self.OGxSeqHashtoKeep[(OGGroup,self.code)] = []
			self.OGxSeqHashProt[(OGGroup,self.code)] = []
		PathtoTemp = self.PathtoTemp

		for record in self.blast:
			if record.descriptions:
				if record.alignments[0].hsps[0].expect < blastCutOff:
					self.getgrpp(record,interestList, PathtoTemp)
				else:
					self.blast.remove(record)
			else:
				self.blast.remove(record)	

	
	def getgrpp(self, record,interestList, PathtoTemp):

		
		if record.descriptions:
			try:
				OGGroup = 'OG5_' +  record.alignments[0].hit_def.split('_')[-1].strip()
		
			except:
				OGGroup = record.alignments[0].hit_def
			
			if OGGroup in interestList:
				
				#out.write(record.query + ',' + OGGroup + '\n')
				self.makeFasta(record.query,OGGroup, PathtoTemp)


	def makeFasta(self, query, OGGroup, PathtoTemp):	#also populates OGxSeqHash
		outname = OGGroup + '_' + self.code  + '_out.txt'
		self.OGlist.append(OGGroup)
		out = open(PathtoTemp + '/' + outname,'a')
		for seq in self.seqs:
			if seq.name == query:
				#print 'adding ' + seq.seq + ' to hash' #NEW 
				if self.is_nuc() == 'nuc':
					self.OGxSeqHash[(OGGroup,self.code)].append(seq)
				else:
					self.OGxSeqHashProt[(OGGroup,self.code)].append(seq)					
				out.write('>' +  seq.name + '\n' + str(seq.seq)  + '\n')
		
		out.close()


###################################################################################################	

###################################################################################################	
#	*Taxon2* getbestblastallseqs, checkFrame, findFS, reversecomp, trans_new, and isOneFrame				  #
#	Translation scripts!																		  #
###################################################################################################
	def getbestblastallseqs(self,table):  				
		for OG in self.OGList:
			self.removeidentical(OG)	
			#print self.OGxSeqHash
			for seq in self.OGxSeqHash[(OG,self.code)]:
				#print seq.seq
				for blast in self.blast:
					if seq.id == blast.query:
						seqOG = 'OG5_' +  blast.alignments[0].title.split('_')[-1]
						self.checkFrame(seq.seq, seq.id, seqOG, blast,table) #checks frame and does the translation
	
	def getbestblastallseqsnt(self,table):  				
		for OG in self.OGList:
			self.removeidentical(OG)	
			#print self.OGxSeqHash
			for seq in self.OGxSeqHashtoKeep[(OG,self.code)]:
				#print seq.seq
				for blast in self.blast:
					if seq.id == blast.query:
						seqOG = 'OG5_' +  blast.alignments[0].title.split('_')[-1]
						self.checkFrame(seq.seq, seq.id, seqOG, blast,table) #checks frame and does the translation	
			os.system('mv ' + self.PathtoOutput + '/fasta2keep/' + OG + '_' + self.code + 'fastatokeep.fas ' + self.PathtoOutput + '/fasta2keep/' + OG + '_' + self.code + '_ntfastatokeep.fas')
			os.system('mv ' + self.PathtoOutput + '/fasta2keep/' + OG + '_' + self.code + 'Transfastatokeep.fas ' + self.PathtoOutput + '/fasta2keep/' + OG + '_' + self.code + 'fastatokeep.fas')
			#copies the translated fasta onto nuc fasta
	
	def checkFrame(self, seq, id, seqOG, blast,table):	

		frame = self.isOneFrame(blast) #this is where it checks to see if it is in one frame
		#print frame
		if frame == 'YY': #multiple reverse reads
			self.reversecomp( seq,id,seqOG,table)
	
		else:	
			
			#try:
				if frame != None:
					if type(frame) == int: #frame is a number, so only one frame, just translate
						self.trans_new( seq, id, table,frame,'all',seqOG)
					else: #frame is a dict, so multiple frames
						self.findFS(frame,table,  seq, id, seqOG)
				else:
					print 'problem is frame = None'
			
	def findFS(self,frameDict,table,  seq, id, seqOG): #deletes one base at frame shift and repeats
		keyList=[]	
		for key in frameDict:
			keyList.append(key)
		newseq = seq[keyList[0]:keyList[1]-1] + seq[keyList[1]:]
		outfile = open(self.PathtoTemp + '/query','w')
		outfile.write('>' + id + '\n')
		outfile.write(str(newseq))
		outfile.close()
		
		db = seqOG
		os.system('formatdb -i ' + self.PathtoTemp + '/'  + seqOG +'  -p T -o F')
		cline = ('blastx -db ' + self.PathtoTemp + '/' + seqOG + ' -query ' + self.PathtoTemp + '/query -num_descriptions 5 -num_alignments 5 -evalue .00001 -num_threads 2 -outfmt 5  -out ' + self.PathtoTemp + '/blast_out.txt')
		os.system(cline)
		try:
			blast_record = NCBIXML.read(open(self.PathtoTemp + '/blast_out.txt'))
		
			if blast_record.descriptions:
				self.checkFrame(newseq, id, seqOG, blast_record,table) 
		except:
			logfile = open(self.PathtoOutput + '/logfile','a')
			logfile.write('Taxon: problem with ' + seqOG + '\n')
			logfile.close()




	def reversecomp(self, seq,id,seqOG,table):

		newseq = seq.reverse_complement()
		outfile = open(self.PathtoTemp + '/query','w')
		outfile.write('>' + id + '\n')
		outfile.write(str(newseq))
		outfile.close()

		db = seqOG
		os.system('formatdb -i ' + self.PathtoTemp + '/'  + seqOG +'  -p T -o F')
		cline = ('blastx -db ' + self.PathtoTemp + '/' + seqOG + ' -query ' + self.PathtoTemp + '/query -num_descriptions 5 -num_alignments 5 -evalue .00001 -num_threads 2 -outfmt 5  -out ' + self.PathtoTemp + '/blast_out2.txt')
		os.system(cline)
		blast_record = NCBIXML.read(open(self.PathtoTemp + '/blast_out2.txt'))
		if blast_record.descriptions:
			self.checkFrame(newseq, id, seqOG, blast_record,table) 
		
	def trans_new(self,   seq, id,table, frame, seqrange, seqOG):
		#print seqOG
		c_uncinata_table = CodonTable(forward_table={
		    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
		    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
		    'TAT': 'Y', 'TAC': 'Y',             'TAG': 'Q',
		    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
    		'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
		    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
		    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
		    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
		    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
		    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
		    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
		    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
		    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
		    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
		start_codons = [ 'ATG'],
		stop_codons = ['TAA' ])


		
		outFileName = self.PathtoOutput + '/' + seqOG + '_' + self.code + 'Trans.fas'	

		outFile=open(outFileName,"a")
	 
		seq_len = len(seq) 
		if seqrange == 'all':
			if frame < 0:
				nuc = seq.reverse_complement()
			else:
				nuc = seq
			if table == 'c_uncinata_table':
				trans = str(nuc[abs(frame)-1:].translate(table=c_uncinata_table)) 
			else:
				trans = str(nuc[abs(frame)-1:].translate(table)) 
	
		splitTrans = trans.split('*')
		ORF1 = max(splitTrans, key = len)
		try:
			ORF2 = re.sub('J','X',ORF1)
			ORF3 = re.sub('B','X',ORF2)
			ORF = re.sub('Z','X',ORF3)
		except:
			print 'problem translating!!'
			ORF = max(splitTrans, key = len)
		outFile.write('>' + id +  '_trans'  + '\n')
		outFile.write(ORF)
		outFile.write('\n')		
		outFile.close()	
		newSeq = SeqRecord(Seq(ORF), id = id)
		#print newSeq.id, newSeq.seq
		self.OGxSeqHashProt[(seqOG,self.code)].append(newSeq)
		

	def trans_newnt(self,   seq, id,table, frame, seqrange, seqOG):
		#print seqOG
		c_uncinata_table = CodonTable(forward_table={
		    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
		    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
		    'TAT': 'Y', 'TAC': 'Y',             'TAG': 'Q',
		    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
    		'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
		    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
		    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
		    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
		    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
		    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
		    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
		    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
		    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
		    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
		start_codons = [ 'ATG'],
		stop_codons = ['TAA' ])


		
		outFileName = self.PathtoOutput + '/fasta2keep/' + seqOG + '_' + self.code + 'Transfastatokeep.fas'	

		outFile=open(outFileName,"a")
	 
		seq_len = len(seq) 
		if seqrange == 'all':
			if frame < 0:
				nuc = seq.reverse_complement()
			else:
				nuc = seq
			if table == 'c_uncinata_table':
				trans = str(nuc[abs(frame)-1:].translate(table=c_uncinata_table)) 
			else:
				trans = str(nuc[abs(frame)-1:].translate(table)) 
	
		splitTrans = trans.split('*')
		ORF1 = max(splitTrans, key = len)
		try:
			ORF2 = re.sub('J','X',ORF1)
			ORF3 = re.sub('B','X',ORF2)
			ORF = re.sub('Z','X',ORF3)
		except:
			print 'problem translating!!'
			ORF = max(splitTrans, key = len)
		outFile.write('>' + id +  '_trans'  + '\n')
		outFile.write(ORF)
		outFile.write('\n')		
		outFile.close()	
		newSeq = SeqRecord(Seq(ORF), id = id)
		#print newSeq.id, newSeq.seq
		self.OGxSeqHashProt[(seqOG,self.code)].append(newSeq)

	def isOneFrame(self,blast_record): #this is where it checks to see if it is in one frame
		framelist=[]
		frameDict={}
		elist = []
		eDict = {}
	

		for hsp in blast_record.alignments[0].hsps:
			frameDict[hsp.query_start] = hsp.frame[0]
			framelist.append(hsp.frame[0])
		framelist.sort()	

		if framelist[0] == framelist[-1]:
			return framelist[0]
		else:
			if framelist[0] > 0: #multiple forward frames
				return frameDict

			else:  
				if framelist[-1] > 0: #two hsps in different directions--take lowest evalue
					for hsp in blast_record.alignments[0].hsps:
						eDict[hsp.expect] = hsp
						elist.append(hsp.expect)
					elist.sort()
					return eDict[elist[0]].frame[0]
				else:#multiple reverse frames
					return "YY"


###################################################################################################	
#	*Taxon3* pairwiseclustalList, removeidentical, sortfasbylen, comparescores and filterlen 			  #
#	Ridding!																					  #
###################################################################################################
	def pairwiseclustalListnt(self): #for yt method, removing alleles before translating
		os.system('mkdir ' + self.PathtoOutput + '/Needlescores')
		os.system('mkdir ' + self.PathtoOutput + '/fasta2keep')
		os.system('mkdir ' + self.PathtoOutput + '/AllSeqFiles')	
	
		for OG in self.OGList:
			file = OG + '_' + self.code
			self.removeidentical(OG)	
			toRemove = []
			toKeep = []
			written = []
			
			fasList = self.sortfasbylennt(OG)	#list of proteins sorted by length and short ones thrown out
			print fasList
			if fasList != []:
				if len(fasList) > 0:
					toKeep.append(fasList.pop(0)) #longest sequence will be kept
					toKeep = self.comparescoresnt(fasList,toKeep,toRemove,file,[0])
					#os.system('mv ' + self.PathtoOutput + '/' + file + ' ' + self.PathtoOutput + '/AllSeqFiles')
					outfas = open(self.PathtoOutput + '/fasta2keep/' + file + 'fastatokeep.fas','a')
					for item in toKeep:
						#print item.seq
						if item.id not in written:
							outfas.write('>' + item.id + '\n' + str(item.seq) + '\n')
							written.append(item.id)
							self.OGxSeqHashtoKeep[(OG,self.code)].append(item)
							os.system('mv ' +  self.PathtoOutput + '/' + file + '_reduced ' + self.PathtoOutput + '/AllSeqFiles')
				else:
					for item in fasList:
						#print item
						self.OGxSeqHashtoKeep[(OG,self.code)].append(item)
					outfasname = file + 'fastatokeep.fas'
					outfas = open(self.PathtoOutput + '/fasta2keep/' + outfasname,'a')
					for seq in self.OGxSeqHashtoKeep[(OG,self.code)]:
						outfas.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
					#os.system('cp '  + self.PathtoOutput + '/' + file + '*reduced ' + self.PathtoOutput + '/fasta2keep/' + outfasname)
					os.system('mv ' +  self.PathtoOutput + '/' + file + '_reduced ' + self.PathtoOutput + '/AllSeqFiles')
			else:
				outfasname = file + 'fastatokeep.fas'
				print 'check that ' + outfasname + ' is nuc.'
				for item in self.OGxSeqHashProt[(OG,self.code)]:
					self.OGxSeqHashtoKeep[(OG,self.code)].append(item)
				outfasname = file + 'fastatokeep.fas'
				outfas = open(self.PathtoOutput + '/fasta2keep/' + outfasname,'a')
				for seq in self.OGxSeqHashtoKeep[(OG,self.code)]:
					outfas.write('>' + seq.id + '\n' + str(seq.seq) + '\n')				
				os.system('mv ' +  self.PathtoOutput + '/' + file + '_reduced ' + self.PathtoOutput + '/AllSeqFiles')
			#print self.OGxSeqHashtoKeep
			outfas.close()
			
			try:
				if os.path.getsize(self.PathtoOutput + '/fasta2keep/' + outfasname) == 0:
					os.system('rm ' + self.PathtoOutput + '/fasta2keep/' + outfasname)
			except:
				flag = 0
			
		os.system('rm ' + self.PathtoOutput + '/OG5_*')
	
	
	def pairwiseclustalList(self):
		os.system('mkdir ' + self.PathtoOutput + '/Needlescores')
		os.system('mkdir ' + self.PathtoOutput + '/fasta2keep')
		os.system('mkdir ' + self.PathtoOutput + '/AllSeqFiles')
		
		
		for OG in self.OGList:
			filename = OG + '_' + self.code
			#self.removeidentical(OG)

			toRemove = []
			toKeep = []
			written = []
	
			fasList = self.sortfasbylen(OG)	#list of proteins sorted by length and short ones thrown out
			if fasList != []:
				if len(fasList) > 0:
					toKeep.append(fasList.pop(0)) #longest sequence will be kept
					toKeep = self.comparescores(fasList,toKeep,toRemove,filename)
					
					#os.system('mv ' + self.PathtoOutput + '/' + file + ' ' + self.PathtoOutput + '/AllSeqFiles')
					
					outfas = open(self.PathtoOutput + '/fasta2keep/' + filename + 'fastatokeep.fas','a')
					for item in toKeep:
					
						#print item.seq
						if item.id not in written:
							outfas.write('>' + item.id + '\n' + str(item.seq) + '\n')
							written.append(item.id)
							self.OGxSeqHashtoKeep[(OG,self.code)].append(item)
							os.system('mv ' +  self.PathtoOutput + '/' + filename + '_reduced ' + self.PathtoOutput + '/AllSeqFiles')
				else:
					outfas = open(self.PathtoOutput + '/fasta2keep/' + filename + 'fastatokeep.fas','a')
					for item in fasList:
						outfas.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
						self.OGxSeqHashtoKeep[(OG,self.code)].append(item)
					#outfasname = fasfile + 'fastatokeep.fas'
					#outfas = open(self.PathtoOutput + '/fasta2keep/' + outfasname,'a')
					#for seq in self.OGxSeqHashtoKeep[(OG,self.code)]:
					#	outfas.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
					#os.system('mv ' +  self.PathtoOutput + '/' + file + '_reduced ' + self.PathtoOutput + '/AllSeqFiles')
				outfas.close()
			#else:
			#	
			#	outfasname = file + 'fastatokeep.fas'
			#	print 'check that ' + outfasname + ' is nuc.'
			#	for item in self.OGxSeqHashProt[(OG,self.code)]:
			#		self.OGxSeqHashtoKeep[(OG,self.code)].append(item)
			#	outfasname = file + 'fastatokeep.fas'
			#	outfas = open(self.PathtoOutput + '/fasta2keep/' + outfasname,'a')
			#	for seq in self.OGxSeqHashtoKeep[(OG,self.code)]:
			#		outfas.write('>' + seq.id + '\n' + str(seq.seq) + '\n')				
			#	os.system('mv ' +  self.PathtoOutput + '/' + file + '_reduced ' + self.PathtoOutput + '/AllSeqFiles')
			##print self.OGxSeqHashtoKeep
			
			
			try:
				if os.path.getsize(self.PathtoOutput + '/fasta2keep/' + outfasname) == 0:
					os.system('rm ' + self.PathtoOutput + '/fasta2keep/' + outfasname)
			except:
				flag = 0
			
		os.system('rm ' + self.PathtoOutput + '/OG5_*')
	def removeidenticalnt(self,OG):
		logfile = open(self.PathtoTemp + '/reductionlog.txt','w')
		logfile.write(OG + ': \n')


		seqlist = []
		outfile = open(self.PathtoOutput + '/' + OG + '_' + self.code + '_reduced','a')
		for seq in self.OGxSeqHash[(OG,self.code)]:
			taxon = seq.id.split('_')[0] + '_' + seq.id.split('_')[1]	

		
			if [taxon,str(seq.seq)] not in seqlist:#if sequence hasn't been seen, add it			
				seqlist.append([taxon,str(seq.seq)])
				outfile.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
			else:
				logfile.write('removed: ' +  seq.id + '\n')
				self.OGxSeqHash[(OG,self.code)].remove(seq)
		#os.system('rm ' + self.PathtoOutput + '/' + f )
		#os.system('mv ' + self.PathtoOutput + '/' + f + '_reduced  '  + self.PathtoOutput + '/' + f)
	
	def removeidentical(self,OG):
		logfile = open(self.PathtoTemp + '/reductionlog.txt','w')
		logfile.write(OG + ': \n')


		seqlist = []
		outfile = open(self.PathtoOutput + '/' + OG + '_' + self.code + '_reduced','a')
		for seq in self.OGxSeqHashProt[(OG,self.code)]:
			taxon = seq.id.split('_')[0] + '_' + seq.id.split('_')[1]	

		
			if [taxon,str(seq.seq)] not in seqlist:#if sequence hasn't been seen, add it			
				seqlist.append([taxon,str(seq.seq)])
				outfile.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
			else:
				logfile.write('removed: ' +  seq.id + '\n')
				self.OGxSeqHashProt[(OG,self.code)].remove(seq)
		#os.system('rm ' + self.PathtoOutput + '/' + f )
		#os.system('mv ' + self.PathtoOutput + '/' + f + '_reduced  '  + self.PathtoOutput + '/' + f)
		
		
	def sortfasbylen(self, OG):
		inFile1 = []
		inFile2 = []
	
		for seq in self.OGxSeqHashProt[(OG,self.code)]:
			#print seq.id
			if str(seq.seq) not in inFile2:		#seqs we have already seen	
				inFile1.append(seq) #list of sequence objects
				inFile2.append(str(seq.seq)) #list of sequences we have seen
		#print inFile1
		#print inFile2
		#outfile=open(self.PathtoTemp + '/sorted.fasta','w') #changed from 'a' to 'w'...should  be better?
	
		for i in range(0, len(inFile1) - 1):
			for j in range(0, len(inFile1) - i - 1):
				if len(inFile1[j].seq) < len(inFile1[j + 1]):
					inFile1[j], inFile1[j + 1] = inFile1[j + 1], inFile1[j]  # sort by length
		if inFile1 != []:
			#print inFile1
			inFile3 = self.filterlen(inFile1)  #get rid of short sequences
	
			#for fas in inFile3:
			#	outfile.write('>' + fas.id + '\n' + str(fas.seq) + '\n') #sortedfasta = file containing sorted sequences
		#print inFile3[0].id
		else:
			inFile3 = inFile1
			if inFile3 == []:
				print 'no Seqs for ' + OG
		return inFile3 #sorted list of seqs
	
	def comparescores(self, fasList,toKeep,toRemove,filename):
		flag = 0
		toRemoveIndex = []
		if len(fasList) > 0:
			for i in range(len(fasList)):				
				out = open(self.PathtoTemp + '/seq1.fasta','w')
				out2 = open(self.PathtoTemp + '/seq2.fasta','w')
				out.write('>' + toKeep[-1].id  + '\n' + str(toKeep[-1].seq) + '\n')
				out2.write('>' + fasList[i].id  + '\n' + str(fasList[i].seq)+ '\n') # compare one to all in the file
				out.close()
				out2.close()
				print filename
				cline = 'EMBOSS-6.5.7/emboss/needle -outfile=' + self.PathtoOutput + '/Needlescores/' + filename + 'needle.txt -asequence=' + self.PathtoTemp + '/seq1.fasta -bsequence=' + self.PathtoTemp + '/seq2.fasta -gapopen=10 -gapextend=0.5' 
				os.system(cline)
						
				
				
				
				x = self.test(filename)

				if x == True:
					flag = 1
				else:	
					toRemoveIndex.append(i)	
			print toRemoveIndex
			toRemoveIndex.reverse()
			print toRemoveIndex
			for index in toRemoveIndex:
				fasList.remove(fasList[index])
				
			if len(fasList) > 0:
				toKeep.append(fasList.pop(0))	
			self.comparescores(fasList,toKeep,toRemove,filename)

		else:
			print 'done'

		settoKeep = set(toKeep)
		toKeep = list(settoKeep)
		settoRemove = set(toRemove)
		toRemove = list(settoRemove)

		for item in toKeep:
			if item in toRemove:
				toKeep.remove(item)
		return toKeep

	def filterlen(self,inFile1):	
		#print inFile1[0].id
		inFile3 = []
		maxlen = len(inFile1[0].seq)
		inFile3.append(inFile1.pop(0))
		
		#print maxlen
		if len(inFile1) > 0:
			for seq in inFile1:
				if len(seq) > 0.7*maxlen:
					inFile3.append(seq)
		return inFile3


	def test(self,f):
		outfinal = open(self.PathtoOutput + '/Needlescores/' + f + 'needlescores.txt','a')
		outscore = open(self.PathtoOutput + '/Needlescores/' + f + 'pairwise_out_scores.csv','a')
		from Bio import AlignIO
		import re
		flag = 1
		infile = open(self.PathtoOutput + '/Needlescores/' + f + 'needle.txt','r')
		midline = ""
		gap = 0.0
		ident = 3.0
		sim = 3.0
		diff = 0.0
	
		for line in infile:
			outfinal.write(line)
			if re.search('# 1: ', line):
				seq1 = str(line.split(':')[1].strip())
			if re.search('# 2: ', line):
				seq2 = str(line.split(':')[1].strip())		
			if re.search('# Length:',line):
				tlengthlist = re.findall('\d+', line)
				tlength = float(tlengthlist[0])
			
		#print line[0]
			if line.strip() == '#=======================================':
				flag = 0
		
			if line.strip() == '#---------------------------------------':
				flag = 1
			
			if flag == 0:
				if line[0] == ' ':
					midline=midline + line.strip()
		outfinal.write('\n')	
		midline = midline + '*'
		
		midline2 = ""
			
		id = 0 
		gp = 0
	
	
	
	
		for char in midline:
			if id < 4:	
				if char == '|':
					id = id + 1
				else: 
					id = 0
			else:
				midline2 = midline2 + char
	
		for char in midline2:
	
			#if gp < 10:	
			if char == ' ':
				gp = gp + 1
				gap = gap + 1
			else:
				gp = 0
			if char == '|':
				ident = ident + 1
				sim = sim + 1
			if char == ':':
				sim = sim + 1
			if char == '.':
				diff = diff + 1
			if char == '*':
				gp = 10
		len = gap + sim + diff
		
		infile.close()
		outscore.write(seq1 + ',' + seq2 + ', len: ' + str(len) + ', tlength: ' + str(tlength) + ', identity: ' + str(ident) + ', similarity: ' + str(sim) + ', len/tlength: ' + str(float(len)/float(tlength)) + ', sim/len: ' + str(float(sim)/float(len)) + '\n')
		##############
		# *Taxon_P1*
		##############

		if float(len)/float(tlength) > 0.4: #length of overlap > 2/3 - too stringent? try 1/4
			if float(sim)/float(len) > 0.70: #num identical or similar at least 3/4 of overlap (similar only means something for aa)
				if float(sim)/float(len) < 0.98: #num identical or similar at most .98 of overlap (rids too similar seqs)
					print 'True'
					return True
