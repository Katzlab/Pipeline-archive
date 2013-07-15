import sys


print '#############################################################################'
print 'Useage is "python renamefasta.py <yourFastaFile>'
print '#############################################################################'
print 'output is renamed file which should be checked before using it in the pipeline'
print '#############################################################################'
print 'This script takes a fasta file downloaded from genbank and renames the sequences with a standard naming format that you enter.\n'
print 'The standard is to rename each sequences with a 2 letter major clade code, a 2 letter minor clade code and a 4 letter species code plus the gi number unique to that sequence.\n'
print 'For example, if you give it the code "Op_me_hsap" for Opisthokont_metazoa_H.sapiens" it will preceed each sequence name with that code, \n'
print 'and will extract the gi number from the standard genbank naming system.  This will only work with fasta files downloaded from genbank, but can be adapted. \n'
print '#############################################################################'

def rename(arg, code):
	out = open(arg + '_renamed.fas','a')
	infile = open(arg,'r')
	for line in infile:
		if line[0] == '>':
			ginum = line.split('|')[1]
			newline = code + ginum + '\n'
		else:
			newline = line

		out.write(newline)


def main():
	
	code = raw_input('What would you like to preceed the unique sequence identifier? i.e. "Op_me_Hsap_" ')
	for arg in sys.argv[1:]:		
		rename(arg, code)

main()
