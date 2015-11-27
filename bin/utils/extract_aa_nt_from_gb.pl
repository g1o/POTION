#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;


chomp $ARGV[0] if defined;
chomp $ARGV[1] if defined;

if ((!defined $ARGV[1])||(!-d $ARGV[0])||(!-d $ARGV[1])||($ARGV[0] eq "-h")||(($ARGV[0] eq "--help"))) { #checking if users are asking for help
  print "\nperl extract_aa_nt_from_gb.pl <directory with genbank files> <directory to write output>\n\n";
  print "The script will read gb files, extract CDS features and translations and write CDS and protein data in two fasta files\n\n";
  print "A third file, linking the common temporary ID used in both files to the old protein/CDS ids is also created\n\n";
  exit();
}


my $dir_path = shift;
my $output_dir = shift;

opendir(DIR,$dir_path);
my @files = readdir(DIR);

foreach my $file (@files) {
  if (($file eq ".") || ($file eq "..")) { next; }
  if ($file !~ /.gb$/) { next; }
  
  print ("$file\n");
  my $in = new Bio::SeqIO(-format => 'genbank',
                          -file => "$dir_path/$file");

  open (OUT, ">", "$output_dir/$file\_proteins.fasta");
  open (OUT_2, ">", "$output_dir/$file\_nt.fasta");
  open (LOG, ">", "$output_dir/$file\_IDs.log");

  my $source = $file;
  $source =~ s/\.gb.*$//g;
  $source =~ s/\s+/_/g;
  $source =~ s/\(/_/g;
  $source =~ s/\)/_/g;
    
  my $i = 1;
  while (my $seq = $in->next_seq()) {
    for my $feat_object ($seq->get_all_SeqFeatures) {

      if ($feat_object->primary_tag ne 'CDS') { next; }
      if ($feat_object->has_tag('pseudo')) { next; }
      
      my $cds_object = $feat_object->spliced_seq;
      my $id = $cds_object->display_id();
      my $cds = $cds_object->seq;
      
      print OUT (">$source\_$i\n");
      print OUT_2 (">$source\_$i\n$cds\n");
      
      if ($feat_object->has_tag('protein_id')) {
        for my $value ($feat_object->get_tag_values('protein_id')) {
          print LOG ("$source\_$i\t$value\n");
          $i++;
        }
      }
      if ($feat_object->has_tag('translation')) {
        for my $value ($feat_object->get_tag_values('translation')) {
          print OUT ("$value\n");
        }
      }
    }
  }
}
