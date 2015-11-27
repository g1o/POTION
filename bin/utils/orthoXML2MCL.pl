#!/usr/bin/env perl


#######################################################################
#######################################################################
# Copyright 2010 Embrapa Informática Agropecuária
# Author: Adhemar Zerlotini Neto
# orthoXML2MCL.pl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
#
# GenerateChromosomesReport.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details <http://www.gnu.org/licenses/>.
#
#######################################################################
#######################################################################


use strict;
use XML::Simple qw(:strict);
use Data::Dumper;


my $xmlFile = shift or die('File name required!!!');

my $in = XMLin($xmlFile, ForceArray => 1, KeyAttr => [] );

#print Dumper($in);


#print $in->{'species'}[0]{'database'}[0]{'genes'}[0]{'gene'}[0]{'protId'};


my $ids = {};
my $flagGene = 0;
my $flagProt = 0;
my @geneIdCount;
my @protIdCount;
foreach my $species (@{$in->{'species'}}) {
	#print $species->{'name'} . "\t" . $species->{'NCBITaxId'} . "\n";
	my $species_name = $species->{'name'}; 
	if($species->{'NCBITaxId'}) {
		$species_name = $species->{'NCBITaxId'};
	} else {
		$species_name =~ s/\s/_/g; 
	}
	foreach my $database (@{$species->{'database'}}) {
		#print $database->{'version'} . "\t" . $database->{'name'} . "\t" . $database->{'protLink'} . "\n";
		foreach my $genes (@{$database->{'genes'}}) {
			foreach my $gene (@{$genes->{'gene'}}) {
				#print $gene->{'id'} . "\t" . $gene->{'geneId'} . "\t" . $gene->{'protId'} . "\n";
				$ids->{$gene->{'id'}}{'geneId'} = $gene->{'geneId'};
				$ids->{$gene->{'id'}}{'protId'} = $gene->{'protId'};
				$ids->{$gene->{'id'}}{'species'} = $species_name;
				
				if( grep( /$gene->{'geneId'}/, @geneIdCount ) ) {
					$flagGene = 1;
				} else {
					#print $gene->{'geneId'} . "\n";
					push @geneIdCount,$gene->{'geneId'};
				}
				if( grep( /$gene->{'protId'}/, @protIdCount ) ) {
					$flagProt = 1;
				} else {
					#print $gene->{'geneId'} . "\n";
					push @protIdCount,$gene->{'protId'};
				} 
				
			}
		}
	}
}

#print Dumper $ids;

my $i = 0;
foreach my $groups (@{$in->{'groups'}}) {
	foreach my $orthologGroup (@{$groups->{'orthologGroup'}}) {
		my $orthomclLine = '';
		if($orthologGroup->{'id'}) {
			$orthomclLine = $orthologGroup->{'id'};
		} else {
			$orthomclLine = "ORTH" . $i;
		}
		my @genesCount;
		my @taxaCount;
		my $auxLine = '';		
		foreach my $geneRef (@{$orthologGroup->{'geneRef'}}) {
			#print $ids->{$geneRef->{'id'}}{'species'};
			unless( grep( /$ids->{$geneRef->{'id'}}{'species'}/, @taxaCount ) ) {
				push @taxaCount,$ids->{$geneRef->{'id'}}{'species'};
			}
			push @genesCount,$geneRef->{'id'};
			if($ids->{$geneRef->{'id'}}{'geneId'} && ! $flagGene) {
				$auxLine .= " " . $ids->{$geneRef->{'id'}}{'geneId'} . "(" . $ids->{$geneRef->{'id'}}{'species'} . ")";
			} elsif($ids->{$geneRef->{'id'}}{'protId'} && ! $flagProt) { 
				$auxLine .= " " . $ids->{$geneRef->{'id'}}{'protId'} . "(" . $ids->{$geneRef->{'id'}}{'species'} . ")";
			} else {
				$auxLine .= " " . $geneRef->{'id'} . "(" . $ids->{$geneRef->{'id'}}{'species'} . ")";
			}
		}
		$orthomclLine .= "(" . int(@genesCount) . " genes," . int(@taxaCount) . " taxa)" . $auxLine . "\n";
		print $orthomclLine;	
		$i++;
	}
}


