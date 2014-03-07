#!/usr/bin/perl -w
use lib "/bioseq/Guidance/";

use strict;
use Storable;

use Guidance;

my $stored_data_file=shift;
my $stored_form_file=shift;
my $Cutoff=shift;

my $vars_ref = retrieve($stored_data_file);
my %VARS = %$vars_ref;

my $form_ref = retrieve($stored_form_file);
my %FORM = %$form_ref;

open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!");

#remove sites with SP-score < Col sp_cutoff
############################################
$VARS{Seq_File_without_low_SP_SEQ}.=".$Cutoff";
$VARS{removed_low_SP_SEQ}.=".$Cutoff";
$VARS{Seq_File_without_low_SP_SEQ_with_Names}=$VARS{Seq_File_without_low_SP_SEQ}.".With_Names";
$VARS{removed_low_SP_SEQ_With_Names}=$VARS{removed_low_SP_SEQ}.".With_Names";


print LOG "Guidance::removeLowSPseq (\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr\",\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}\",$Cutoff,\"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}\");\n";
my @ans=Guidance::removeLowSPseq ("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr","$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}",$Cutoff,"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}");
print "ANS:",join("",@ans),"\n";
if (-s  "$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}" > 0) # NOT EMPTY
{
	my @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ_with_Names}");
	if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ_with_Names}\") failed:",join("",@ans),"\n");}
}

if (-s "$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}" > 0) # Seq were removed
{
	my  @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{removed_low_SP_SEQ_With_Names}");
	if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ_With_Names}\") failed:".join("",@ans)."\n");}
}


# Update the output page
#######################################
open (OUTPUT,"$VARS{WorkingDir}$VARS{output_page}");
my @out=<OUTPUT>;
close (OUTPUT);
open (OUTPUT,">$VARS{WorkingDir}$VARS{output_page}");
my $Remove_Seq_Section=0;
foreach my $line (@out)
{
    if ($line=~/Remove unreliable sequences below confidence score/)
	{
		$Remove_Seq_Section=1;
		print OUTPUT $line;
	}
	elsif (($line=~/form/) and $Remove_Seq_Section==1 and ($line!~/form.data/))
	{
		print OUTPUT $line;
		if ($FORM{'Redirect_From_MAFFT'}==1) {print_message_to_output("<A HREF='$VARS{Seq_File_without_low_SP_SEQ_with_Names}' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below $Cutoff)</A><font size=-1> (see list of removed sequences <A HREF='$VARS{removed_low_SP_SEQ_With_Names}' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"var answer = confirm('ATTENTION: Running GUIDANCE on the confidently-aligned sequences only, ignores the parameters used for the original run on MAFFT server. It is therefore recommended to adjust these parameters or aligning the confidently-aligned sequences on MAFFT server and run GUIDANCE again from there');if (answer){window.open('http://guidance.tau.ac.il/index_rerun.php?run=$VARS{run_number}&file=$VARS{Seq_File_without_low_SP_SEQ_with_Names}')}\"><br>");}
		else {print_message_to_output("<A HREF='$VARS{Seq_File_without_low_SP_SEQ_with_Names}' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below $Cutoff)</A><font size=-1> (see list of removed sequences <A HREF='$VARS{removed_low_SP_SEQ_With_Names}' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"window.open('http://guidance.tau.ac.il/index_rerun.php?run=$VARS{run_number}&file=$VARS{Seq_File_without_low_SP_SEQ_with_Names}')\"><br>");}
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
