IVA Publication
===============

This repository contains all supplementary data, scripts and
intermediate data necessary to reproduce the results of
the IVA manuscript (submitted).

For all of the scripts in this repository to work, you will
need IVA (and its dependencies) installed, and also
the Python package [Fastaq] [fastaq] installed.
This README
assumes that the reader has read the publication and supplementary
material.

------------------------------------------------------------------------------

Description of directories
==========================

Assemblies
----------
The `Assemblies/` directory contains two directories:
`HIV/` and `Flu/`, each of which contains one
directory per sample. Each sample directory has a
directory for each of the assemblers IVA, IVA.c5r2,
PRICE, Trinity, Inchworm, VICUNA.80 and VICUNA.90.

The files within each assembly directory are:

  * `contigs.fasta` - FASTA file of contigs made by the assembler
  * `iva_qc.contig_placement.pdf` - contig layout diagram made by `iva_qc`
  * `iva_qc.reads_mapped_to_assembly.bam.flagstat` - file of mapping
    statistics onto the assembly. Made by running `samtools flagstat` as part of `iva_qc`.
  * `iva_qc.reads_mapped_to_ref.bam.flagstat` - file of mapping
    statistics onto the reference genome. Made by running `samtools flagstat` as part of `iva_qc`.
  * `iva_qc.stats.txt` - file of QC metrics output by the script `iva_qc`
  * `seed_reads.fa` - FASTA file of seed reads used as input to PRICE (only
     applicable to PRICE directory).



Data
----
The files in the `Data/` directory are described in the supplementary material.



Figure 1
--------
The `Figure_1` directory contains data and the script needed to make
Figure 1 of the paper. Run `./make_figure_1.py` to regenerate
the file `figure_1.svg`.

The script `bam_to_cov_and_snps.pl` is also included. This was run on
the BAM file output by `iva_qc` to make the files `contigs.iva.fa.read_depth.gz`
and `contigs.iva.fa.read_depth.gz.tbi`, which are needed by
`make_figure_1.py`.



Scripts
-------
The `Scripts/` directory contains custom analysis scripts used
for the publication.

To make all plots (except Figure 1, which is described above)
and regenerate the QC summary files, run

    ./Scripts/gather_results.py Assemblies/ plot
    mv plot.flu.tsv Data/table.S5.qc_summary.flu.tsv
    mv plot.hiv.tsv Data/table.S6.qc_summary.hiv.tsv

To make the plots of CPU and RAM usage, run

    cd Data
    ../Scripts/make_resources_plots.R table.S7.resources_flu.tsv table.S8.resources_hiv.tsv resources

This script outputs the rows of the outliers that were removed from the
box plots of wall clock and total CPU time, then makes the plots.

The usage of the remaining scripts is described in the supplementary material.


  [fastaq]: https://github.com/sanger-pathogens/Fastaq
