#!/usr/bin/perl -w
use lib "/bioseq/Guidance/";

use strict;
use Storable;

use Guidance;

my $stored_data_file=shift;

my $Cutoff=shift;

my $vars_ref = retrieve($stored_data_file);
my %VARS = %$vars_ref;

open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!");

#remove sites with SP-score < Col sp_cutoff
############################################
$VARS{Alignment_File_without_low_SP_Col}.=".$Cutoff";
$VARS{removed_low_SP_SITE}.=".$Cutoff";

print LOG "Guidance::removeLowSPsites (\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr\",\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}\",$Cutoff,\"$VARS{WorkingDir}$VARS{removed_low_SP_SITE}\");\n";
my @ans=Guidance::removeLowSPsites ("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr","$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}",$Cutoff,"$VARS{WorkingDir}$VARS{removed_low_SP_SITE}");
if ($ans[0]eq "OK")
  {
    $VARS{REMOVED_SITES}=$ans[1];
    $VARS{MSA_LENGTH}=$ans[2];
  }
print LOG "REMOVED_SITES:$VARS{REMOVED_SITES}\n";
print LOG "MSA_LENGTH:$VARS{MSA_LENGTH}\n";
$VARS{Alignment_File_without_low_SP_Col_with_Names}=$VARS{Alignment_File_without_low_SP_Col}.".With_Names";
if (-s "$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}" > 0) # Not EMPTY
{
	my @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col_with_Names}");
	if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col_with_Names}\") failed:",join("",@ans),"\n");}
      }
# Update the output page
#######################################
open (OUTPUT,"$VARS{WorkingDir}$VARS{output_page}");
my @out=<OUTPUT>;
close (OUTPUT);
open (OUTPUT,">$VARS{WorkingDir}$VARS{output_page}");
my $Remove_Pos_Section=0;
foreach my $line (@out)
  {
    if ($line=~/Remove unreliable columns below confidence score/)
      {
	$Remove_Pos_Section=1;
	print OUTPUT $line;
      }
    elsif (($line=~/form/) and $Remove_Pos_Section==1)
      {
	print OUTPUT $line;
	print_message_to_output("<A HREF='$VARS{Alignment_File_without_low_SP_Col_with_Names}' TARGET=_blank>The MSA after the removal of unreliable columns (below $Cutoff)</A><font size=-1> (see list of removed columns <A HREF='$VARS{removed_low_SP_SITE}' TARGET=_blank>here</A>)</font><br>"); 
	$Remove_Pos_Section=0;
      }
    else
      {
	print OUTPUT $line;
      }
  }
close (OUTPUT);


#---------------------------------------------
sub print_message_to_output{
#---------------------------------------------
    my $msg = shift;
    print OUTPUT "\n<ul><li>$msg</li></ul>\n";
}
