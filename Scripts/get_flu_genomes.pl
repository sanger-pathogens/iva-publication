#!/usr/bin/env perl
use strict;
use warnings;
use File::Slurp;
use List::Util qw(shuffle);

system("wget ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/genomeset.dat") and die $!;
my $lines = read_file('genomeset.dat');
my @genome_blocks = split(/\n\n/, $lines);

my @flu_A_genomes;
my @flu_B_genomes;

for my $block (@genome_blocks) {
    my @records = split /\n/, $block;
    next if (scalar @records) != 8;
    my $ok = 1;
    my %segments;
    my @ids;
    my %a_or_b;
    for my $i (0..$#records) {
        $records[$i] = [split /\t/, $records[$i]];
        if ($records[$i][1] ne 'Human') {
            $ok = 0;
            last;
        }
        if ($records[$i][7] =~ /Influenza ([AB]) virus/) {
            $a_or_b{$1} = 1;
        }
        else {
            $ok = 0;
            last;
        }
        $segments{$records[$i][2]} = 1;
        push @ids, $records[$i][0];
    }

    next if not $ok or (scalar keys %a_or_b != 1) or (not defined $a_or_b{"A"} and not defined $a_or_b{"B"});

    for my $i (1..8) {
        unless (exists $segments{$i}) {
            $ok = 0;
            last;
        }
    }

    next if not $ok;

    if (defined $a_or_b{"A"}) {
        push @flu_A_genomes, join(' ', @ids);
    }
    elsif (defined $a_or_b{"B"}) {
        push @flu_B_genomes, join(' ', @ids);
    }


}

@flu_A_genomes = shuffle(@flu_A_genomes);
@flu_B_genomes = shuffle(@flu_B_genomes);

print join("\n", @flu_A_genomes[0..99]), "\n";
print join("\n", @flu_B_genomes[0..99]), "\n";
