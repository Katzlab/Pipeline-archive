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
Author: Jessica Grant.  Email questions/bug reports to jgrant@smith.edu
Please cite: (TBD) 