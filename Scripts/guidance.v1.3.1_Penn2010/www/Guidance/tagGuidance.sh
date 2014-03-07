#!/bin/csh
if ($#argv == 0) then
	echo "USAGE: tagGuidance.sh <version number>"
else
	set v = $argv[1];
	mkdir tags/guidance.v$v;
	svn add tags/guidance.v$v;
	svn cp trunk/Makefile tags/guidance.v$v/;
	mkdir tags/guidance.v$v/libs;
	svn add tags/guidance.v$v/libs;
	svn cp trunk/libs/Makefile tags/guidance.v$v/libs/;
	svn cp trunk/libs/phylogeny tags/guidance.v$v/libs/;
	mkdir tags/guidance.v$v/programs;
	svn add tags/guidance.v$v/programs;
	svn cp trunk/programs/Makefile tags/guidance.v$v/programs/;
	svn cp trunk/programs/Makefile.generic tags/guidance.v$v/programs/;
	svn cp trunk/programs/semphy tags/guidance.v$v/programs/;
	svn cp trunk/programs/msa_set_score tags/guidance.v$v/programs/;
	svn cp trunk/programs/isEqualTree tags/guidance.v$v/programs/;
	mkdir tags/guidance.v$v/www;
	svn add tags/guidance.v$v/www;
	svn cp trunk/www/bioSequence_scripts_and_constants tags/guidance.v$v/www/;
	svn cp trunk/www/Guidance tags/guidance.v$v/www/;
	svn cp trunk/www/Selecton tags/guidance.v$v/www/;
	svn mv tags/guidance.v$v/www/Guidance/README tags/guidance.v$v/;
endif

# Don't forget to remove irrelevant libs and programs from the makefiles 
