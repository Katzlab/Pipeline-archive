Pipeline
========
unzip the package, pipelinepkg.zip, and you will find two folders:

Scripts: contains the python scripts needed to run the pipeline, including some of the stand-alone software packages the pipeline relies upon.
Files: contains the data files needed to run the pipeline.  It also has the test run file "test"  To run your own pipeline, a similar file listing your OGs of interest will need to be placed in the Files folder.

The Files folder is currently set up to run a test case, the phylogenetic analysis of Monosiga ovata.  For a different analysis, remove the Monosiga files from the folders TaxonDataFiles and BlastFiles and replace with your own Taxon and Blast files (The Blast files can be made from within the script, but it may be faster (or at least less user input during the pipeline analysis) to run the stand-along blast from ncbi (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) separately.

Due to size limitations on github, we are unable to provide a necessary file - allOG5Files - but we will be happy to send this to anyone who wants it.  Or, these data can be downloaded from OrthoMCL (orthomcl.org).

All scripts and folders must remain in the paths they come in, or else the scripts will need to be reworked.
The following dependencies must be installed on your computer:

python - (up to version 2.7.5 - Python 3 is not supported)
biopython (http://biopython.org/wiki/Download) 
dendropy (http://pythonhosted.org/DendroPy/)
RAxML-PTHREADS (http://sco.h-its.org/exelixis/software.html) called as 'raxml'.  To call differently, change the following line in Gene/__init__.py to reflect your raxml call:
       raxCL = ('raxml -T 4 -s ' + file + ' -n ' + self.OG + '_outrax.tree -f d -m PROTGAMMALG') 
mafft (http://mafft.cbrc.jp/alignment/software/)
formatdb and blast ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Other software called (i.e. Guidance, needle) are in the Scripts package and will be called with that path.  

Please cite:
Guidance: Penn O, Privman E, Ashkenazy H, Landan G, Graur D, Pupko T: GUIDANCE: a web server for assessing alignment confidence scores. Nucleic Acids Res 2010, 38:W23-W28.
Needle: Rice,P. Longden,I. and Bleasby,A. "EMBOSS: The European Molecular Biology Open Software Suite" Trends in Genetics June 2000, vol 16, No 6. pp.276-277
========
To run:

cd into the Scripts folder and run phylopipe.py as follows:

python phylopipe.py

You will be asked to ether the name of your OGofInterestFile, the text file with a list of OGs (see test for an example)
This OGofInterestFile must be in the Files folder.
========
Output Files:

The output will be found in a folder with the same name as your OGofInterest file.  In it there will be a zip file of all output from all stages of the pipeline.
There will also be three uncompressed folders with the most important output:

RenamedAlignments - full alignments with full names
BestRaxTreesRenamed - single gene trees made from the alignments in RenamedAlignments
AlignmentsforConcatenation - alignments with paralogs removed and named for concatenation
========
Test case:
The test case is set up to run the pipeline on one taxon - Monosiga ovata - and two genes. 

The taxon data is in Files/TaxonDataFiles
The blast data is in Files/BlastFiles (preblasted against all data downloaded from orthoMCL - this will be created in the script if not found)
The data downloaded from orthoMCL is in two files in Files/allOG5Files
the text file test is a list of the OGs to be analyzed - the same ones as are in the allOG5Files folder.

cd into Scripts and type <python phylopipe.py>

you will be asked for your OGofInterestList.  Type in <test.txt> (the name of the file in Files with a list of OGs.)

The script looks into the TaxonDataFiles to find the taxa to be added.  In the test case, there is only one taxon, Monosiga ovata.

A new folder is made named for the output, named for the OGofInterestList, test.txt

The pipeline first collects the orthologs, refines the data using pairwise alignements, then runs guidance to create alignments for each gene.

The intermediate results are in directories within test.txt/Output/ which will be compressed when the pipeline completes successfully

Output files will be in the test.txt/Output/

The test has been run and the results are in the folder 'test_results' in the github repository
 
========
Author: Jessica Grant.  Email questions/bug reports to jgrant@smith.edu
Please cite: (TBD) 