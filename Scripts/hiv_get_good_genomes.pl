#!/usr/bin/env perl
use strict;
use warnings;
use LWP::Simple;

@ARGV == 1 or die "usage: $0 <HIV1_COM_2012_genome_DNA.fasta>

Finds genomes with all 9 expected HIV genes.
Prints list of IDs to stdout
";


my @ids = @{alignemnt_file_to_ids($ARGV[0])};

for my $id (@ids) {
    my $fname = "tmp.$id.gb";
    download_gb($id, $fname);
    if (gb_has_all_genes($fname)) {
        print "$id\n";
    }
    unlink $fname;
}


sub alignemnt_file_to_ids {
    my $filename = shift;
    my @ids;
    open F, $filename or die $!;
    while (my $line = <F>) {
        next unless $line =~ /^>/;
        next if $line =~ /^>CPZ/;
        chomp $line;
        my @a = split(/\./, $line);
        push(@ids, $a[-1]);
    }
    close F;
    return \@ids;
}


sub looks_like_genbank {
    my $filename = shift;
    open F, $filename or die $!;
    my $line = <F>;
    close F;
    return $line =~/^LOCUS/;
}


sub download_gb {
    my $id = shift;
    my $outfile = shift;
    my $maxtries = 5;
    my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=text&id=$id";
    foreach (1..$maxtries) {
        getstore($url, $outfile);    
        return if looks_like_genbank($outfile);
        sleep(3);
    }
    die "Error downloading $id";
}


sub gb_has_all_genes {
    my $filename = shift;
    my %wanted_genes = (
        env => 0,
        gag => 0,
        nef => 0,
        pol => 0,
        rev => 0,
        tat => 0,
        vif => 0,
        vpr => 0,
        vpu => 0,
    );

    open F, $filename or die $!;

    while (<F>) {
        chomp;
        if (/^                     \/gene="(.+)"$/) {
            delete $wanted_genes{$1} if exists $wanted_genes{$1};
        } 
    }

    close F or die $!;
    return 0 == scalar keys %wanted_genes ? 1 : 0;
}

