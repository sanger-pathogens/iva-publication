#!/usr/bin/env python3

import sys
import os
import copy
import math
import argparse
import fastaq
import iva


assemblers = ['IVA', 'IVA.c5r2', 'PRICE', 'Inchworm', 'Trinity', 'VICUNA.80', 'VICUNA.90']
assemblers_main = ['IVA', 'PRICE', 'Trinity', 'VICUNA.80']
assemblers_short = {
    'IVA': 'IVA',
    'PRICE': 'PRI',
    'VICUNA.80': 'VIC',
    'Trinity': 'Tri'
}

qc_stats = [
    'trimmed_reads',
    'reads_map_ref',
    'reads_map_assembly',
    'ref_bases',
    'ref_sequences',
    'ref_bases_assembled',
    'ref_sequences_assembled',
    'ref_sequences_assembled_ok',
    'ref_bases_assembler_missed',
    'assembly_bases',
    'assembly_bases_in_ref',
    'assembly_contigs',
    'assembly_contigs_hit_ref',
    'assembly_bases_reads_disagree',
    'assembly_sum_longest_match_each_segment',
    'cds_number',
    'cds_assembled',
    'cds_assembled_ok',
    'gage_Missing_Reference_Bases',
    'gage_Missing_Assembly_Bases',
    'gage_Missing_Assembly_Contigs',
    'gage_Duplicated_Reference_Bases',
    'gage_Compressed_Reference_Bases',
    'gage_Bad_Trim',
    'gage_Avg_Idy',
    'gage_SNPs',
    'gage_Indels_<_5bp',
    'gage_Indels_>=_5',
    'gage_Inversions',
    'gage_Relocation',
    'gage_Translocation',
    'ratt_elements_found',
    'ratt_elements_transferred',
    'ratt_elements_transferred_partially',
    'ratt_elements_split',
    'ratt_parts_of_elements_not_transferred',
    'ratt_elements_not_transferred',
    'ratt_gene_models_to_transfer',
    'ratt_gene_models_transferred',
    'ratt_gene_models_transferred_partially',
    'ratt_exons_not_transferred_from_partial_matches',
    'ratt_gene_models_not_transferred',
]


class Assembly:
    def __init__(self, sample_id, assembler, contigs_file, qc_file, flagstat_assembly, flagstat_ref):
        self.id =  sample_id
        self.assembler = assembler
        self.contigs_file = contigs_file
        self.flagstat_assembly = flagstat_assembly
        self.flagstat_ref = flagstat_ref
        self.qc_stats_file = qc_file
        self.qc = {x: -1 for x in qc_stats}
        self._load_qc_stats()
        self._load_map_stats()

    def _load_qc_stats(self):
        f = fastaq.utils.open_file_read(self.qc_stats_file)
        for line in f:
            key, value = line.rstrip().split('\t')
            if key.startswith('ref_EMBL'):
                continue
            assert key in self.qc
            if '.' in value:
                self.qc[key] = float(value)
            elif value.isdigit():
                self.qc[key] = int(value)

        fastaq.utils.close(f)

    def _mapstats_to_stats(self, filename):
        f = fastaq.utils.open_file_read(filename)
        line = f.readline()
        if ' in total ' in line:
            total_reads = int(line.split()[0])
        else:
            print("Didn't get total reads from ", filename)
            print('line:', line)
            sys.exit(1)

        line = f.readline()
        line = f.readline()
        
        if ' mapped ' in line:
            mapped_reads = int(line.split()[0])
        else:
            print("Didn't get mapped reads from ", filename)
            print('line:', line)
            sys.exit(1)

        fastaq.utils.close(f)
        return total_reads, mapped_reads

    def _load_map_stats(self):
        total_assembly, mapped_assembly = self._mapstats_to_stats(self.flagstat_assembly) 
        total_ref, mapped_ref = self._mapstats_to_stats(self.flagstat_ref)
        assert total_assembly == total_ref
        self.qc['trimmed_reads'] = total_ref
        self.qc['reads_map_ref'] = mapped_ref
        self.qc['reads_map_assembly'] = mapped_assembly

    def to_tsv_data_line(self):
        return '\t'.join([
            self.id,
            self.assembler,
            '\t'.join([str(self.qc[x]) for x in qc_stats])
        ])


class Assemblies:
    def __init__(self, root_dir, outprefix):
        self.root_dir = os.path.abspath(root_dir)
        self.outprefix = outprefix
        self.tsv_data_file = self.outprefix + '.tsv'
        self.r_script = self.outprefix + '.R'
        self.tex_stats = {}
        self.samples = 0
        self._get_assembly_data()
        self._get_assemblers()
        self._filter_assemblies()

    def _get_assembly_data(self):
        self.assemblies = {}
        for sample in os.listdir(self.root_dir):
            sample_dir = os.path.join(self.root_dir, sample)
            assert sample not in self.assemblies
            self.assemblies[sample] = []
            self.samples += 1

            for assembler in assemblers:
                contigs_file = os.path.join(sample_dir, assembler, 'contigs.fasta')
                qc_file = os.path.join(sample_dir, assembler, 'iva_qc.stats.txt')
                flagstat_assembly = os.path.join(sample_dir, assembler, 'iva_qc.reads_mapped_to_assembly.bam.flagstat')
                flagstat_ref = os.path.join(sample_dir, assembler, 'iva_qc.reads_mapped_to_ref.bam.flagstat')
                ok = True
                for filename in [contigs_file, qc_file, flagstat_ref, flagstat_assembly]:
                    if not os.path.exists(filename):
                        print('Warning, missing file(s) for', sample_dir, assembler, filename, file=sys.stderr)
                        ok = False

                if ok:
                    a = Assembly(sample, assembler, contigs_file, qc_file, flagstat_assembly, flagstat_ref)
                    self.assemblies[sample].append(a)


    def _get_assemblers(self):
        all_assembler_tuples = set()
        for assembly_id in self.assemblies:
            assemblers = tuple(sorted([x.assembler for x in self.assemblies[assembly_id]]))
            all_assembler_tuples.add(assemblers)

        biggest_tuple_size = max([len(x) for x in all_assembler_tuples])
        biggest_tuples = [x for x in all_assembler_tuples if len(x) == biggest_tuple_size]
        number_of_biggest = len(biggest_tuples)
        if number_of_biggest != 1:
            print('Error! conflicting assembler names. Cannot continue', file=sys.stderr)

        self.assemblers = set(biggest_tuples[0])


    def _filter_assemblies(self):
        to_delete = set()
        for assembly_id in self.assemblies:
            assemblers = set([x.assembler for x in self.assemblies[assembly_id]])
            if assemblers != self.assemblers:
                print('Ignoring ID', assembly_id, 'because not all assemblers present. Needed:', self.assemblers, ' got:', assemblers)
                to_delete.add(assembly_id)

        for assembly_id in to_delete:
            del self.assemblies[assembly_id]

        print('Using', len(self.assemblies), 'of', self.samples, 'samples')


    def _write_tsv_data_file(self):
        f = fastaq.utils.open_file_write(self.tsv_data_file)
        print('ENA_ID', 'Assembler', '\t'.join(qc_stats), sep='\t', file=f)
        for assembly_id in self.assemblies:
            for assembly in self.assemblies[assembly_id]:
                print(assembly.to_tsv_data_line(), file=f)
        fastaq.utils.close(f)


    def _get_stats_for_one_assembler(self, assembler, stat):
        stats = []
        for sample in sorted(self.assemblies):
            for a in self.assemblies[sample]:
                if a.assembler == assembler:
                    stats.append(a.qc[stat])
        return stats


    def _get_stats_for_assemblers(self, stat):
        return {assembler: self._get_stats_for_one_assembler(assembler, stat) for assembler in self.assemblers}

    def _dict_to_means(self, d, dp=1):
        return {key: round(sum(d[key]) / len(d[key]), dp) for key in d}

    def _dict_to_stdevs(self, d):
        return {key: self._list_to_stdev(d[key]) for key in d}

    def _list_to_stdev(self, l):
        mean = sum(l) / len(l)
        variance = sum([(x-mean)*(x-mean) for x in l]) / len(l)
        return round(math.sqrt(variance), 2)

    def _list_to_median(self, l):
        a = copy.copy(l)
        a.sort()
        return round(a[int(len(a)/2)], 2)

    def _dict_to_medians(self, d):
        return {key: self._list_to_median(d[key]) for key in d}

    def _gather_tex_stats(self):
        self.tex_stats['samples'] = self.samples
        ref_seqs_assembled_ok = self._get_stats_for_assemblers('ref_sequences_assembled_ok')
        ref_sequences = self._get_stats_for_assemblers('ref_sequences')
        for k in ref_sequences:
            assert len(set(ref_sequences[k])) == 1
        ref_sequences = ref_sequences['IVA'][0]
        good_assemblies = {k: round(100 * len([1 for x in ref_seqs_assembled_ok[k] if x == ref_sequences]) / self.samples, 1) for k in self.assemblers}

        ref_seqs_assembled = self._get_stats_for_assemblers('ref_sequences_assembled')
        ref_bases = self._get_stats_for_assemblers('ref_bases')
        ref_bases_assembled = self._get_stats_for_assemblers('ref_bases_assembled')
        ref_bases_assembler_missed = self._get_stats_for_assemblers('ref_bases_assembler_missed')

        ratt_elements_transferred = self._get_stats_for_assemblers('ratt_elements_transferred')
        ratt_elements_found = self._get_stats_for_assemblers('ratt_elements_found')
        ratt_elements_trans_perc = {k: [100 * ratt_elements_transferred[k][i] / ratt_elements_found[k][i] for i in range(len(ratt_elements_found[k]))] for k in self.assemblers}

        ratt_gene_models_to_transfer = self._get_stats_for_assemblers('ratt_gene_models_to_transfer')
        ratt_gene_models_transferred = self._get_stats_for_assemblers('ratt_gene_models_transferred')
        ratt_genes_trans_perc = {k: [100 * ratt_gene_models_transferred[k][i] / ratt_gene_models_to_transfer[k][i] for  i in range(len(ratt_gene_models_to_transfer[k]))] for k in self.assemblers}

        percent_assembled = {}
        for a in self.assemblers:
            percent_assembled[a] = [100 * ref_bases_assembled[a][i] / ref_bases[a][i] for i in range(len(ref_bases[a]))]


        gage_dups = self._get_stats_for_assemblers('gage_Duplicated_Reference_Bases')
        gage_dup_rate = {}
        for a in self.assemblers:
            gage_dup_rate[a] = [gage_dups[a][i] / ref_bases[a][i] for i in range(len(ref_bases[a]))]

        gage_inv = self._get_stats_for_assemblers('gage_Inversions')
        gage_rel = self._get_stats_for_assemblers('gage_Relocation')
        gage_tran = self._get_stats_for_assemblers('gage_Translocation')
        gage_errors = {}
        for k in self.assemblers:
            gage_errors[k] = [gage_inv[k][i] + gage_rel[k][i] + gage_tran[k][i] for i in range(len(gage_inv[k]))]
        total_bases_assembled = {k: sum(ref_bases_assembled[k]) for k in self.assemblers}
        gage_errors_per_kb_assembled = {k: round(1000 * sum(gage_errors[k]) / total_bases_assembled[k], 3) for k in self.assemblers}
        gage_total_errors = {k: sum(gage_errors[k]) for k in gage_errors}

        gage_avg_idy = self._get_stats_for_assemblers('gage_Avg_Idy')

        assembly_sum_longest_match_each_segment = self._get_stats_for_assemblers('assembly_sum_longest_match_each_segment')
        assembly_sum_longest_match_each_segment_perc = {}
        assembly_bases_in_ref = self._get_stats_for_assemblers('assembly_bases_in_ref')
        assembly_bases_in_ref_perc = {}
        for a in self.assemblers:
            assembly_sum_longest_match_each_segment_perc[a] = [100 * assembly_sum_longest_match_each_segment[a][i] / ref_bases[a][i] for i in range(len(ref_bases[a]))]
            assembly_bases_in_ref_perc[a] = [100 * assembly_bases_in_ref[a][i] / ref_bases[a][i] for i in range(len(ref_bases[a]))]
 
        

        self.tex_stats['gage_Avg_Idy'] = self._dict_to_means(gage_avg_idy)
        self.tex_stats['gage_Avg_Idy_stdev'] = self._dict_to_stdevs(gage_avg_idy)
        self.tex_stats['percent_bases_assembled'] = self._dict_to_means(percent_assembled)
        self.tex_stats['percent_bases_assembled_stdev'] = self._dict_to_stdevs(percent_assembled)
        self.tex_stats['assembly_sum_longest_match_each_segment_perc'] = self._dict_to_means(assembly_sum_longest_match_each_segment_perc)
        self.tex_stats['assembly_sum_longest_match_each_segment_perc_stdev'] = self._dict_to_stdevs(assembly_sum_longest_match_each_segment_perc)
        self.tex_stats['assembly_bases_in_ref_perc'] = self._dict_to_means(assembly_bases_in_ref_perc)
        self.tex_stats['assembly_bases_in_ref_perc_stdev'] = self._dict_to_stdevs(assembly_bases_in_ref_perc)
        self.tex_stats['ratt_elements_trans_perc'] = self._dict_to_means(ratt_elements_trans_perc)
        self.tex_stats['ratt_elements_trans_perc_stdev'] = self._dict_to_stdevs(ratt_elements_trans_perc)
        self.tex_stats['ratt_genes_trans_perc'] = self._dict_to_means(ratt_genes_trans_perc)
        self.tex_stats['ratt_genes_trans_perc_stdev'] = self._dict_to_stdevs(ratt_genes_trans_perc)
        self.tex_stats['good_assemblies'] = good_assemblies
        self.tex_stats['gage_dup_rate'] = self._dict_to_means(gage_dup_rate)
        self.tex_stats['gage_dup_rate_stdev'] = self._dict_to_stdevs(gage_dup_rate)
        self.tex_stats['gage_errors'] = self._dict_to_means(gage_errors)
        self.tex_stats['gage_total_errors'] = gage_total_errors
        self.tex_stats['gage_errors_per_kb_assembled'] = gage_errors_per_kb_assembled

    def _write_R_script(self):
        f = fastaq.utils.open_file_write(self.r_script)
        print('library(ggplot2)',
              'd = read.csv(file="' + self.tsv_data_file + r'''", sep="\t", header=T)''',
              'pre="' + self.outprefix + '."',
              'x_labels = c(' + ', '.join(['"' + a + '"' for a in assemblers]) + ')',
              r'''
ggplot() + geom_boxplot(data=d, mapping=aes(x=Assembler, y=100 * ref_bases_assembled / ref_bases)) +
    xlab("Assembler") +
    scale_x_discrete(labels=x_labels, limits=x_labels) +
    ylab("Per cent of reference assembled") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename=paste(pre, "suppl.percent_of_ref_assembled.pdf", sep=""), scale=0.55)


ggplot() + geom_boxplot(data=d, mapping=aes(x=Assembler, y=ref_sequences_assembled_ok)) +
    xlab("Assembler") +
    ylab("Uniquely assembled segments") +
    scale_x_discrete(labels=x_labels, limits=x_labels) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename=paste(pre,"suppl.uniquely_assembled_segments.pdf", sep=""), scale=0.55)
''', 
               '\n'.join([a + '_ok = dim(d[(d$Assembler == "' + a + '") & (d$ref_sequences_assembled_ok == 1),])[1]' for a in assemblers]),
               'total = dim(d[(d$Assembler == "IVA"), ])[1]',
               'assemblers = c("' + '", "'.join(assemblers) + '")',
               'counts = c(' + ', '.join([a + '_ok' for a in assemblers]) + r''')
df = data.frame(assemblers, counts)
ggplot() + geom_bar(data=df, stat="identity", mapping=aes(x=assemblers, y=100*counts/total), fill="white", colour="black") +
    xlab("Assembler") +
    scale_x_discrete(labels=x_labels, limits=x_labels) +
    ylab("Assembled into a unique contig (%)") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename=paste(pre, "suppl.genome_one_unique_contig.pdf", sep=""), scale=0.55)


ggplot() + geom_boxplot(data=d, mapping=aes(x=Assembler,y=100*ratt_elements_transferred/ratt_elements_found)) +
    xlab("Assembler") +
    scale_x_discrete(labels=x_labels, limits=x_labels) +
    ylab("Annotation features transferred (%)") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename=paste(pre, "suppl.ratt_elements_transferred.pdf", sep=""), scale=0.55)
''', sep='\n', file=f)
        fastaq.utils.close(f)

    def run(self):
        self._write_tsv_data_file()
        self._write_R_script()
        iva.common.syscall('R CMD BATCH ' + self.r_script)
        self._gather_tex_stats()


parser = argparse.ArgumentParser(
    usage = '%(prog)s <rootdir> <outprefix>')
parser.add_argument('rootdir', help='Name of root Assemblies/ directory')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()



hiv = Assemblies(os.path.join(options.rootdir, 'HIV'), options.outprefix + '.hiv')
hiv.run()
flu = Assemblies(os.path.join(options.rootdir, 'Flu'), options.outprefix + '.flu')
flu.run()


def tex_stat_line(stat, stat_str, a):
    fields = [stat_str] + [a.tex_stats[stat][assembler] for assembler in assemblers]
    return ' & '.join([str(x) for x in fields]) + r''' \\'''
         
         
def tex_stat_line_stdev(stat, stat_str, a):
    fields = [stat_str] + [str(a.tex_stats[stat][assembler]) + ' (' + str(a.tex_stats[stat + '_stdev'][assembler]) + ')' for assembler in assemblers]
    return ' & '.join(fields) + r''' \\'''


def write_tex_table(filename, assemblers, assemblies):
    f = fastaq.utils.open_file_write(filename)
    print(r'''\begin{tabular}{l|''' + 'r'*len(assemblers) + r'''}\toprule''', file=f)
    print(' & ' + ' & '.join(assemblers) + r''' \\ \midrule''', file=f)
    print(tex_stat_line('good_assemblies', r'''Ideal assemblies (\%)$^1$''', assemblies), file=f)
    print(tex_stat_line_stdev('percent_bases_assembled', r'''Mean reference bases assembled (\%)''', assemblies), file=f)
    print(tex_stat_line_stdev('assembly_sum_longest_match_each_segment_perc', r'''Longest contig(s) sum (\% of reference)''', assemblies), file=f)
    print(tex_stat_line_stdev('assembly_bases_in_ref_perc', r'''Assembly length (\% of reference)''', assemblies), file=f)
    print(tex_stat_line_stdev('gage_dup_rate', r'''Mean duplication rate$^2$''', assemblies), file=f)
    print(tex_stat_line_stdev('ratt_elements_trans_perc', r'''Mean \% annotation transferred''', assemblies), file=f)
    print(tex_stat_line('gage_total_errors', r'''Total assembly errors$^3$''', assemblies), file=f)
    print(tex_stat_line_stdev('gage_Avg_Idy', r'''Mean per-sample identity to reference (\%)$^4$''', assemblies), r''' \bottomrule''', file=f)
    print(r'''\end{tabular}''', file=f)
    fastaq.utils.close(f)


write_tex_table(options.outprefix + '.suppl.table.hiv.tex', assemblers, hiv)
write_tex_table(options.outprefix + '.suppl.table.flu.tex', assemblers, flu)

all_data_tsv = options.outprefix + '.all.tsv'
f = fastaq.utils.open_file_read(options.outprefix + '.hiv.tsv')
hiv_data = f.readlines()
fastaq.utils.close(f)
f = fastaq.utils.open_file_read(options.outprefix + '.flu.tsv')
flu_data = f.readlines()
fastaq.utils.close(f)

f = fastaq.utils.open_file_write(all_data_tsv)
print('Organism', hiv_data[0], sep='\t', end='', file=f)
for line in hiv_data[1:]:
    print('HIV-1', line, sep='\t', end='', file=f)
for line in flu_data[1:]:
    print('Influenza', line, sep='\t', end='', file=f)
fastaq.utils.close(f)

fig2_r_script = options.outprefix + '.fig2.R'
f = fastaq.utils.open_file_write(fig2_r_script)

print('library(ggplot2)',
      'd = read.csv(file="' + all_data_tsv + r'''", sep="\t", header=T)''',
      'd$Organism = factor(d$Organism, levels=c("HIV-1", "Influenza"))',
      'pre="' + options.outprefix + '"',
      'x_labels_main = c(' + ', '.join(['"' + a + '"' for a in assemblers_main]) + ')',
      'x_labels_main_short = c(' + ', '.join(['"' + a + '"' for a in [assemblers_short[x] for x in assemblers_main]]) + ')',
      'x_labels = c(' + ', '.join(['"' + a + '"' for a in assemblers]) + ')',
      r'''
ggplot(data=subset(d[100 * d$assembly_bases_in_ref / d$ref_bases < 750,], Assembler %in% c(''' + ', '.join(['"' + x + '"' for x in assemblers_main]) + r''')), aes(x=Assembler, y=100 * assembly_bases_in_ref / ref_bases)) +
    geom_hline(aes(yintercept=100), linetype="longdash") +
    #geom_boxplot(colour="red", outlier.colour="red") + # was used to color new results in red for reviewers
    geom_boxplot() +
    theme_classic() +
    ylab("Assembly length / reference length (%)") +
    theme(strip.background = element_rect(colour="black", fill="white")) +
    scale_x_discrete(labels=x_labels_main_short, limits=x_labels_main) +
    scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,650)) +
    facet_grid(. ~ Organism)''',
      r'''ggsave(filename=paste(pre,".fig2b.pdf", sep=""), scale=0.5)''',
      sep='\n', file=f)


print(r'''

ggplot(data=d[100 * d$assembly_bases_in_ref / d$ref_bases < 750,], aes(x=Assembler, y=100 * assembly_bases_in_ref / ref_bases)) +
    geom_hline(aes(yintercept=100), linetype="longdash") +
    geom_boxplot() +
    theme_classic() +
    ylab("Assembly length / reference length (%)") +
    theme(strip.background = element_rect(colour="black", fill="white"), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    scale_x_discrete(labels=x_labels, limits=x_labels) +
    scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,650)) +
    facet_grid(. ~ Organism)''',
      r'''ggsave(filename=paste(pre,".suppl.assembly_length.pdf", sep=""), scale=0.55)''',
      sep='\n', file=f)

print(r'''
ggplot(data=subset(d, Assembler %in% c(''' + ', '.join(['"' + x + '"' for x in assemblers_main]) + r''')), aes(x=Assembler, y=100 * assembly_sum_longest_match_each_segment / ref_bases)) +
    geom_hline(aes(yintercept=100), linetype="longdash") +
    #geom_boxplot(colour="red", outlier.colour="red") + # was used to color new results in red for reviewers
    geom_boxplot() +
    theme_classic() +
    ylab("Longest contig(s) / reference length (%)") +
    theme(strip.background = element_rect(colour="black", fill="white")) +
    scale_x_discrete(labels=x_labels_main_short, limits=x_labels_main) +
    #scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,650)) +
    facet_grid(. ~ Organism)''',
      r'''ggsave(filename=paste(pre,".fig2a.pdf", sep=""), scale=0.5)''',
      sep='\n', file=f)

print(r'''
ggplot(data=d, aes(x=Assembler, y=100 * assembly_sum_longest_match_each_segment / ref_bases)) +
    geom_hline(aes(yintercept=100), linetype="longdash") +
    geom_boxplot() +
    theme_classic() +
    ylab("Longest contig(s) / reference length (%)") +
    theme(strip.background = element_rect(colour="black", fill="white"), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    scale_x_discrete(labels=x_labels, limits=x_labels) +
    #scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,650)) +
    facet_grid(. ~ Organism)''',
      r'''ggsave(filename=paste(pre,".suppl.longest_contig.pdf", sep=""), scale=0.55)''',
      sep='\n', file=f)


fastaq.utils.close(f)
iva.common.syscall('R CMD BATCH ' + fig2_r_script)
os.unlink('Rplots.pdf')
os.unlink(options.outprefix + '.flu.suppl.genome_one_unique_contig.pdf')
os.unlink(options.outprefix + '.hiv.suppl.uniquely_assembled_segments.pdf')
