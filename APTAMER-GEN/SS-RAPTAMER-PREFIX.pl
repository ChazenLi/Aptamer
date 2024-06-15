#!/usr/bin/perl -w
# Illustrate how to build nucleic acids.
use strict;
use warnings;
use NucleicAcidDiscoveryScript;

my $filename = "E:/APTAMER-GEN/s.txt";
my $saveDir = "E:/APTAMER-GEN/";
my $prefix = "TTT";

open my $fh, '<', $filename or die "Could not open '$filename' $!\n";

while (my $line = <$fh>) {
    chomp $line;
    print "$line\n";
    
    my $document = Mdm::Document::Create();
    
    # Create and append to a DNA duplex molecule.
    my $naType = Mdm::rnaSingleStrand;
    my $nucleicAcidMolecule =
      $document->CreateNucleicAcid($prefix, $naType, Mdm::aHelix );
    my $senseChain1 = $nucleicAcidMolecule->Chains->Item(0);
    
    my $na3Prime    = $senseChain1->NucleotideAt3PrimeTerminus;
    $na3Prime = $document->GrowNucleicAcid($na3Prime, $line, Mdm::aHelix );
    
    
    # Ligate the two molecules using the sense strand.
    $document->DeselectAll();
    $na3Prime->Select;
    
    # Mutate first sense nucleotide to thymidine.
    $document->DeselectAll();
    $senseChain1->NucleotideAt5PrimeTerminus->Select;
    $document->MutateNucleotide( $senseChain1->NucleotideAt5PrimeTerminus,
        Mdm::thymidine );
    
    $document->Save("$saveDir$prefix-$line.pdb","pdb");
}