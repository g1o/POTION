#!/usr/local/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "\nperl split_multi_organism_gb.pl <genbank file containing more than one species>\n\n";

if (!$ARGV[0]) {
  die ($usage);
}

my $infile = $ARGV[0];
my $in = new Bio::SeqIO(-format => 'genbank',
                        -file => $infile);



my $seq;
my %organisms;

while($seq = $in->next_seq()){
  for my $feat_object ($seq->get_SeqFeatures) {
#    print "primary tag: ", $feat_object->primary_tag, "\n";
    for my $tag ($feat_object->get_all_tags) {
      if ($tag eq "organism") {
        for my $specie ($feat_object->get_tag_values($tag)) {
#          print ("$specie\n");
          if (defined $organisms{$specie}) {
#            print ("Opa!\n");
  	  }
	  else {
	    $organisms{$specie} = 1;
	  }
        }
      }
      for my $value ($feat_object->get_tag_values($tag)) {
#         print "    value: ", $value, "\n";
      }
    }
  }
}

foreach my $key (keys %organisms) {
  my $outfile = $key;
  $outfile =~ s/\s+/_/g;
  $outfile =~ s/-/_/g;
  $outfile =~ s/\\/_/g;
  $outfile =~ s/\//_/g;
  $outfile =~ s/\(/_/g;
  $outfile =~ s/\)/_/g;
  print ("$outfile\n");
  my $in = new Bio::SeqIO(-format => 'genbank',
                          -file => $infile);
  my $seqout = Bio::SeqIO->new(-file   => ">$outfile.gb",
                               -format => 'genbank');

  my $seq;
  while ($seq = $in->next_seq()){
    for my $feat_object ($seq->get_SeqFeatures) {
#    print "primary tag: ", $feat_object->primary_tag, "\n";
      for my $tag ($feat_object->get_all_tags) {
        if ($tag eq "organism") {
          for my $specie ($feat_object->get_tag_values($tag)) {
            print ("\t->$specie\n");
            if ($key eq $specie) {
              print ("\t\t!!!$key\n");
              $seqout->write_seq($seq);
            }
          }
        }
      }
    }
  }
} 
