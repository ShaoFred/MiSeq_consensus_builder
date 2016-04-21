#!/usr/bin/perl
################################################
# primerID_ZS_cutoff_consensus.pl
#
# Wei Shao
# Advanced Biomedical Computing Center
# Leidos Biomedical Research, Inc
# Frederick National Laboratory for Cancer Research
# Maryland, USA
#
# 04-20-2016
#
# Usage: primerID_ZS_cutoff_consensus.pl <blast_file>
# 
# Make sure $front_anchor (line 40) and $back_anchor (line 41) match the upstream and downstream fragments of NNNNNNNNNN in the sample reads
# The coverage cutoff model is based on Zhou, S. et. al. 2015, J Virol 89:8540–8555
# 
############################## #######################


use locale;
use strict;
use POSIX qw/ceil/;

my ($input_blast, $output_file, $Dogtag_list, $block, @blastblock);
my ($blastresult, $r21, $r22, $sampleID, $id, $rest, @alignblock, @query_subj_lines, @refseqline, @sampleline);
my ($s, $each_query, @each_compare, @aligned_lines, $line,@single_refline, @single_sampleline,$first_ref,$last_ref,$first_sample);
my ($last_sample, $ref_seq, $sample_seq, @ref_seq, $firstdiff_length, $Dogtag, %dogtag, $base, @base_pos, $p);
my ($highest_dogtag_count, $consens_cutoff,  $read_number, $ID3, $seqgroup_ref);
my (%seqgroup, %Tag_ID_seq_hash_hash, $Tag, $consensus_seq);
my (%dogtag_list, $DogTag_ID_hash_hash_ref, %DogTag_ID_hash_hash);
my ($front_anchor, $back_anchor, $dogtag_cover);


$input_blast = shift;

$output_file = $input_blast."_consensus.fas";
$Dogtag_list = $input_blast."_primerID_list.xls";

$front_anchor = "TCCA";  ## change to match your sample
$back_anchor = "CTAGCA"; ## change to match your sample

$dogtag_cover = 0.8;    ## used in sub DogTagConsensus

my ($sec,$min,$hour,$day,$month,@rest) =   localtime(time);

open (OUT, ">$output_file");

open (IN, "$input_blast") || die "cannot open $input_blast\n";;

while (<IN>) {

    if (/Query=\s+/ .. /Effective\s+search\s+space\s+used:\s+\d+/) {

	$block .=$_;
	    }
    else {
	push @blastblock, $block if ($block);
        $block = "";
        ### Note: the last alignment block will be skipped because nothing after "/Effective/" patten.
        ### if we add any letters after it, then it will be processed. 
    }
}

 foreach $blastresult (@blastblock) {

  next if ($blastresult =~ /^BLASTN/);

   if (($blastresult =~ /(^.*) rank/) || ($blastresult =~ /(^.*)(\s+)(\S+)/) ){
	
        ($r21, $r22) = split(/\+\@M/, $3);
  
        $sampleID = $1.$r21;
        $sampleID =~ s/\s+//g;
	
    }

   ($id, $rest) = split('length=', $sampleID); ## This is 454 pattern
  
    @query_subj_lines = split (/>/, $blastresult);

    @refseqline =();
    @sampleline = ();

    for ($s = 1; $s < scalar@query_subj_lines; $s++) { ## only one query block.
        $each_query = $query_subj_lines[$s]; 
	
        @each_compare = split (/Strand=/,$each_query); 
  	       
	@aligned_lines = split(/\n/, $each_compare[1]);
        shift@aligned_lines;
	    
	foreach $line (@aligned_lines) {

	   if ($line =~ /Query/) {
		
	       $line =~ s/Query\s+//;
               @single_sampleline = split(/\s+/, $line);
	       @sampleline =(@sampleline, @single_sampleline);
	       @single_sampleline=();
		    
	   }
	   if ($line =~ /Sbjct/) {
	       $line =~ s/Sbjct\s+//;

               @single_refline = split(/\s+/, $line);
	       @refseqline = (@refseqline,@single_refline); 
	       @single_refline =();		    
	   }
		
	} ## foreach $line (@aligned_lines)
 
    $first_ref = shift(@refseqline);
			
    $last_ref = pop(@refseqline);
    $first_sample = shift(@sampleline);
    $last_sample = pop(@sampleline);
    
    $ref_seq = "@refseqline";
    $sample_seq = "@sampleline";	

    $ref_seq =~ s/\s+|\d+//g;
    $sample_seq =~ s/\s+|\d+//g;
    $ref_seq = uc($ref_seq);
    $sample_seq = uc($sample_seq);
     
    next if ($last_ref < $first_ref); ## There should be no reversed seq here.	
    next if ($last_ref - $first_ref <200);  ## There are some very short alignment, probalby from phix:  
        
    ## initialization for each sample.
			
    @ref_seq = ();  

    $firstdiff_length = $first_ref-1; # the length of ref fragment that is not aligned.

	push @ref_seq,  $first_ref, $ref_seq, $sample_seq;  ## @ref_seq should be better named as @ref_sample.  ### Jan 18, 2013

        $Dogtag =  ScanDogTag($id, $ref_seq, $sample_seq, $firstdiff_length);
   
	$dogtag{$id} = $Dogtag."%@ref_seq" if ($Dogtag);
	$dogtag_list{$id} = $Dogtag  if ($Dogtag);

     last; 
   } 

} ## foreach $blastresult; An alignment has been processed here. 

$s = 0;
foreach $id (keys %dogtag) {

   ($Dogtag, $rest) = split(/%/, $dogtag{$id} );

if ( ($Dogtag =~ /-/) || ($Dogtag =~ /N/i)) {
	delete $dogtag{$id}; 
	delete $dogtag_list{$id}; 

    }
}

($sec,$min,$hour,$day,$month,@rest) =   localtime(time);
# printf qq{Time:\t%02d:%02d:%02d\n}, $hour, $min, $sec;

$highest_dogtag_count = MostAbundentDogtag(\%dogtag_list); ## 06-18-2015

$consens_cutoff = RS_Cutoff($highest_dogtag_count);  ## 06-18-2015

 $consens_cutoff = 5 if ($consens_cutoff < 5);  ## 06-18-2015

open (Dout, ">$Dogtag_list");   ### moved here 0114-2015
print Dout "cutoff = $consens_cutoff\n";
foreach $id (keys %dogtag_list) {
    print Dout "$id\t$dogtag_list{$id}\n";  ## used for dogtag freq calculation in dogtag_distribution.pl 

}


close(Dout);

$DogTag_ID_hash_hash_ref = DogtagGrouper(\%dogtag);   ## 06-24-2015
%DogTag_ID_hash_hash = %$DogTag_ID_hash_hash_ref;   ## change name to hash_hash 



## This section delete any dogtags shared by <5 sequences and delete sequences that with thouse dogtags. 
foreach $Tag (keys %DogTag_ID_hash_hash ) {
  
     $read_number = keys %{$DogTag_ID_hash_hash{$Tag} };
	if ($read_number <$consens_cutoff ) {   ## 06-18-2015

	    delete $DogTag_ID_hash_hash{$Tag};
        }
}


foreach $Tag (keys %DogTag_ID_hash_hash ) {
   $read_number = keys %{$DogTag_ID_hash_hash{$Tag} };


   $seqgroup_ref = RemoveInsertion(\%{$DogTag_ID_hash_hash{$Tag} } );   
   %seqgroup = %$seqgroup_ref;
   $Tag_ID_seq_hash_hash{$Tag} = {%seqgroup};

    $consensus_seq  = DogTagConsensus ($seqgroup_ref); #### remark this line if only need dogtag distribution, not dogtag consensus. 

   next unless ($consensus_seq);     ### 12-1-2015

     print OUT ">$Tag $read_number\n$consensus_seq\n";



}

##################################################################
# sub ScanDogTag
###################################################################

sub ScanDogTag {
    my ($sampleid, $ref_seq1, $sample_seq1, $length_diff) = @_;

     my ($i, $refbases, $reffront, $refafter,  $samplebases, $samplefront, $sampleafter, $ref_length, $found_tag);
  
    $ref_length = length($ref_seq1);
 
    $refbases = "";
    $samplebases ="";
    $found_tag = 0; 

     for ($i=340; $i < $ref_length - 10; $i++) { 

	    $refbases = substr($ref_seq1, $i, 10);
	    
	    $samplebases = substr($sample_seq1, $i, 10);
    
	    if ($refbases eq "NNNNNNNNNN") {
  
		if ( ( substr($sample_seq1, $i-4, 4) eq $front_anchor) && ($samplebases !~ /-/) && (substr($sample_seq1, $i+10, 6 ) eq $back_anchor ) ){
		    

		    $found_tag = 1;	
	            last;
	        }
            }
	    
      }


    return $samplebases if $found_tag;
}
#######################################################
# sub DogtagGrouper
# This sub finds IDs of sequences that have the same dog tag. 
########################################################
sub DogtagGrouper {

 my  $dog_tag_ref = shift; 
 my (%dog_tag_list, $id_s, $tag, %dogtag_ID_hash_hash, $seqpair);   #### change %dogtag_ID_hash_array to %dogtag_ID_hash_hash

 %dog_tag_list = %$dog_tag_ref;



 foreach $id_s (keys %dog_tag_list) {
   

   ($tag, $seqpair) = split(/%/, $dog_tag_list{$id_s} );
 

   $dogtag_ID_hash_hash{$tag}{$id_s} = $seqpair;

 }


 return \%dogtag_ID_hash_hash;

}


############################################################
# sub RemoveInsertion
# This subroutine removes insertions from sample sequences. 
############################################################

sub RemoveInsertion {
    
    my  $dog_id_ref  = shift;
    my (%dog_seqpair, $id1, $temp_seq, $ref_seq1,$sample_seq1, $number_deletions ); ## 05-15-2013 added $number_deletions to compensate 5' truncations
    my (%id_seq1, $first_ref_pos, $i);
    my (@temp_array);

    %dog_seqpair = %$dog_id_ref;  

    foreach $id1 (keys %dog_seqpair) {
	$temp_seq ="";

        
	@temp_array = split(/\s+/, $dog_seqpair{$id1});
        $first_ref_pos = $temp_array[0];		## Jan 18, 2013
	$ref_seq1 = $temp_array[1];    ## Jan 18, 2013
	$sample_seq1 = $temp_array[2];	 ## Jan 18, 2013
	
	foreach ($i = 0; $i < length$ref_seq1; $i++) {
	      if (substr($ref_seq1, $i, 1)  ne "-") {
                   $temp_seq .= substr($sample_seq1, $i, 1);
		
              }
       }     
       $id_seq1{$id1} = substr($temp_seq, (22-$first_ref_pos));

       $number_deletions = "-" x ($first_ref_pos - 1);   ## 05-15-2013, to compensate 5' trucations so that all sequences have the same length.
       $id_seq1{$id1} = $number_deletions.$temp_seq;     ## 05-15-2013, added $number_deletions to compensate 5' trucations so that all sequences have the same length.

#	print OUT  ">$id1\n$id_seq1{$id1}\n"; 	## print all sequences that have a particular dogtag.
    }
 
 return \%id_seq1;
}


################################################
# sub DogTagConsensus
################################################

sub DogTagConsensus {


  my $hash_ref = shift;
  my %ID_seq = %$hash_ref;

 my (@seq, $seq_read, $base, $A, $C, $G, $T, $minus, $most, $total_base, $seq_length, $most_base, $consensus_seq, $j);
  my $gaps = 0;      ### 12-1-2015

@seq = values (%ID_seq);

  $consensus_seq = "";

$seq_length = length$seq[0];



for ($j = 0; $j<= $seq_length-1; $j++) {
    $A = 0;
    $C = 0;
    $G = 0;
    $T = 0;
    $minus = 0;
    $most = 0;
    $total_base = 0;
    foreach $seq_read (@seq) {
	
	$base = substr($seq_read, $j, 1);
       
	$A++ if ($base eq "A");
        $C++ if ($base eq "C");
	$G++ if ($base eq "G");
	$T++ if ($base eq "T");
	$minus++ if ($base eq "-");
    }
 #   print "$i A $A\t C $C\t G $G\t T $T\t minus $minus\n\n";
    $total_base = $A + $C + $G + $T + $minus;
    $most = $A;
    $most_base = "A";
    if ($C > $most) {
	$most = $C;
        $most_base = "C";
    }
    if ($G > $most) {
	$most = $G;
	$most_base = "G";
    }
    if ($T > $most) {
	$most = $T;
	$most_base = "T";
    }
    if ($minus > $most) {
	$most = $minus;
	$most_base = "-";
    }

    if ($most/$total_base < $dogtag_cover) {
	$most_base = "-";
	$gaps++;      ## 12-1-2015
    }
    elsif ($total_base > $most) {
	
	$most_base = lc$most_base;
    }
	
    $consensus_seq .= $most_base;
   
    
}

  $consensus_seq = "" if ($gaps >1);
  return $consensus_seq if ($consensus_seq);
}

#######################################
# sub MostAbundentDogtag
#####################################

sub MostAbundentDogtag {

    my $dogtaglist1 = shift;

    my (%listDogtags, $dgID, %dogtag_count1, @dg_counts_unsort, @dg_count_sort, $biggest_dogtag, $primerID);

    %listDogtags = %$dogtaglist1;

    foreach $dgID (values %listDogtags) {

          if (exists $dogtag_count1{$dgID} ) {
	      $dogtag_count1{$dgID}++;
	  }
	  else {
	      $dogtag_count1{$dgID} = 1;
          }
    }
    @dg_counts_unsort = values(%dogtag_count1);
    @dg_count_sort = sort {$b <=> $a } @dg_counts_unsort; ## sort descending. 
    $biggest_dogtag = $dg_count_sort[0];
    
    return $biggest_dogtag;

}
    

########################################
# sub RS_Cutoff
#######################################

sub RS_Cutoff { 

my $most = shift;

my ($cutoff);

if ($most <=10) {
    $cutoff = 2;
}

elsif ($most <=8500) {
    $cutoff = -1.24*10**-21*($most**6) + 3.53*10**-17*($most**5) - 3.90*10**-13*($most**4) + 2.12*10**-9*($most**3) - 6.06*10**-6*($most**2) + 1.80*10**-2*$most + 3.15;
}
else {
    $cutoff = 0.0079 *$most + 9.4869;
}

$cutoff = ceil($cutoff);

$cutoff = 2 if ($cutoff <3);

return $cutoff;
    
}
