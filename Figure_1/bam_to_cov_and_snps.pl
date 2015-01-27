#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


@ARGV == 3 or die "usage: $0 <in.bam> <ref.fa> <outprefix>";

my $bam = $ARGV[0];
my $ref = $ARGV[1];
my $out = $ARGV[2];


unless (-e "$ref.fai") {
    system("samtools faidx $ref.fai") and die "error indexing reference";
}


my %ref_lengths;
my %cov;
my %snps;
open F, "$ref.fai" or die $!;
while (<F>) {
    my ($name, $length) = split /\t/;
    $ref_lengths{$name} = $length;
    $cov{$name} = [(0) x $length];
    $snps{$name} = [(0) x $length];
}
close F or die $!;


open F, "samtools mpileup -d 999999999 -B -q 20 -Q 20 -f $ref $bam | cut -f 1,2,4,5 |" or die $!;
while (my $line = <F>) {
    chomp $line;
    my ($name, $pos, $depth, $bases) = split(/\t/, $line);
    $cov{$name}->[$pos-1] = $depth;
    $snps{$name}->[$pos-1] = pileup_str_to_snp_count(\$bases);
}
close F or die $!;


open F, "|bgzip -c > $out.gz" or die $!;
for my $name (keys %ref_lengths) {
    for my $i (0 .. ($ref_lengths{$name}-1)) {
        print F join("\t",
            $name,
            $i + 1,
            $cov{$name}[$i],
            $snps{$name}[$i],
        ), "\n";
    }
}
close F;


system("tabix -b 2 -e 2 $out.gz") and die "Error running tabix";


sub pileup_str_to_snp_count{
    my $str = shift;
    $$str =~ s/\^.//g;
    remove_indels($str, "+");
    remove_indels($str, "-");
    $$str =~ s/[,.\$]//g;
    return length $$str;
}



sub remove_indels {
    my $str = shift;
    my $char = shift;

    while(1) {
        my $i = index($$str,$char);
        last if $i == -1;
        my $length_to_remove;
        if (substr($$str, $i+1) =~ /^(\d+)\D/) {
            $length_to_remove = $1 + 1 + length($1);
        }
        else {
            die "Error remove_indels ", $$str;
        }

        substr($$str, $i, $length_to_remove) = '';
    }
}


