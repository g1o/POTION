use strict;
use warnings;

my $infile = $ARGV[0];

my $version = $ARGV[1];

open (IN, "<$infile");

while (my $line = <IN>) {
  chomp $line;
  my @aux = split (/\s/, $line);
  my $group_id = shift @aux;
  chop $group_id; #removing ":"
  my $total_genes = $#aux + 1;
  my $total_species = 0;
  my %tmp; #to calculate total species
  foreach my $element (@aux) {
    my @aux_2 = split (/\|/, $element);
#    print "$aux_2[0]\t$aux_2[1]\n";
#    my $a = <STDIN>;
    if (defined ($tmp{$aux_2[0]})) {   
    }
    else {
      $total_species++;
      $tmp{$aux_2[0]} = 1;
    }
  }
  print "$group_id($total_genes genes,$total_species taxa):\t";
  for (my $i = 0; $i <= $#aux - 1; $i++) {
    my $element = $aux[$i];
    my @aux_2 = split (/\|/, $element);
    print "$aux_2[1]($aux_2[0]) ";
  }
  my $element = $aux[$#aux];
  my @aux_2 = split (/\|/, $element);
  print "$aux_2[1]($aux_2[0])\n";
}
#LEISH1000: lma|LmjF.34.1640 lma|LmjF.36.1270 lma|LmjF.34.1920 lma|LmjF.34.1760 lma|LmjF.34.1780 lma|LmjF.34.1860 lma|LmjF.34.1940 lma|LmjF.34.1700 lma|LmjF.34.1820 lma|LmjF.34.1980 lma|LmjF.34.1800 lma|LmjF.34.1960 lma|LmjF.34.1880 lma|LmjF.34.1900 lma|LmjF.34.1740 lma|LmjF.34.1660 lma|LmjF.34.1680 lma|LmjF.34.1720 lma|LmjF.34.1620 lma|LmjF.34.0960 lma|LmjF.34.1080 lma|LmjF.34.1600 lma|LmjF.34.1560 lma|LmjF.34.1580 lbr|LbrM.18.0460 lbr|LbrM.20.0780 lbr|LbrM.20.4290 lbr|LbrM.20.4300 lbr|LbrM.20.4310 lbr|LbrM.20.4340 lbr|LbrM.20.4320 lbr|LbrM.20.1070 lbr|LbrM.20.1080 lbr|LbrM.20.1480 lbr|LbrM.20.2390 lbr|LbrM.20.2410 lbr|LbrM.20.2370 lbr|LbrM.20.0800 lbr|LbrM.20.1060 lbr|LbrM.08.0290 lbr|LbrM.08.0310 lbr|LbrM.08.0300 lbr|LbrM.08.1060 lbr|LbrM.08.1050 lbr|LbrM.20.0790 lbr|LbrM.08.1130 lbr|LbrM.08.0670 lbr|LbrM.13.1330 lme|LmxM.33.0960 lme|LmxM.33.0961 lme|LmxM.33.1560 lme|LmxM.33.1920e lme|LmxM.33.1920a lme|LmxM.33.1920d lme|LmxM.33.1920b lme|LmxM.33.1920c lme|LmxM.33.1920 lme|LmxM.33.1980 lme|LmxM.36.1270 lme|LmxM.33.1800 lme|LmxM.33.1820 lme|LmxM.33.1580 lme|LmxM.33.1900 lme|LmxM.33.1740 lme|LmxM.33.1720 lme|LmxM.33.1720b lme|LmxM.33.1900a lme|LmxM.33.1720a lme|LmxM.33.1720c lme|LmxM.33.1721 lme|LmxM.33.1560a lme|LmxM.33.1725 lin|LinJ.34.1710 lta|LtaP34.1810 lin|LinJ.34.1680 lin|LinJ.34.1150 lin|LinJ.34.1700 ldo|LdBPK_341150.1 ldo|LdBPK_341690.1 lin|LinJ.34.1690 ldo|LdBPK_341700.1
#

#open (IN, "<$infile");

#while (my $line = <IN>) {
#  chomp $line;
#  my @aux = split (/\t/, $line);
#  my $groupid = shift @aux;
#  my $tmp = $aux[0];
#  $tmp =~ s/\s+$//g;
#  $tmp =~ s/^\s+//g;
#  my @genes = split(/\s+/, $tmp);
#  print "$groupid\n\n";
#  my $tmp_id = parse_cluster_id($groupid);
#  print "$tmp_id:";
#  foreach my $gene (@genes) {
#    next if ($gene eq "");
#    print "$gene\n";
#    my $a = <STDIN>;
#    my ($geneid, $species) = parse_gene_id($gene);
#    print " $species|$geneid";
#  }
#  print "\n";
#  my $a = <STDIN>;
#}

#close IN;

#sub parse_cluster_id {
#  my $line = shift;
#  $line =~ s/:$//;
#  $line =~ s/\(\s*\d*\s*gene[s]?\s*,\d*\s*tax(a|on)\s*\)$//;
#  return $line;
#}

#sub parse_gene_id {
#  my @aux = split (/\(/, $_[0]);
#  my $specie = $aux[1];
#  $specie =~ s/\)//g;
#  return ($aux[0], $specie);  # aux[0] has the if o the gene
#}

