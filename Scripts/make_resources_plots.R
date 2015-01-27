#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(TRUE)
flu_infile = args[1]
hiv_infile = args[2]
outprefix = args[3]
cpu_cutoff = 100 * 60 * 60
wall_cutoff = 50 * 60 * 60

flu=read.csv(file=flu_infile, sep="\t", header=T)
hiv=read.csv(file=hiv_infile, sep="\t", header=T)

flu$Organism = rep("Influenza", dim(flu)[1])
hiv$Organism = rep("HIV", dim(hiv)[1])

all=rbind(hiv, flu)

skipped_all_cpu = all[all$CPU.s. > cpu_cutoff,]
skipped_all_cpu
all_cpu = all[all$CPU.s. <= cpu_cutoff,]

skipped_all_wall = all[all$Wall.s. > wall_cutoff,]
skipped_all_wall
all_wall = all[all$Wall.s. <= wall_cutoff,]


ggplot(data=all_cpu, aes(x=Assembler, y=CPU.s./(60*60))) +
    geom_boxplot() +
    theme_classic() + 
    theme(strip.background = element_rect(colour="black", fill="white"), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    ylab("Total CPU time (hours)") +
    facet_grid(. ~ Organism)
ggsave(filename=paste(outprefix, ".cpu.pdf", sep=""), scale=0.45)

ggplot(data=all_wall, aes(x=Assembler, y=Wall.s./(60*60))) +
    geom_boxplot() +
    theme_classic() + 
    theme(strip.background = element_rect(colour="black", fill="white"), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    ylab("Wall clock time (hours)") +
    facet_grid(. ~ Organism)
ggsave(filename=paste(outprefix, ".wall.pdf", sep=""), scale=0.45)


ggplot(data=all, aes(x=Assembler, y=RAM.GB.)) +
    geom_boxplot() +
    theme_classic() + 
    theme(strip.background = element_rect(colour="black", fill="white"), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    ylab("Peak RAM (GB)") +
    facet_grid(. ~ Organism)
ggsave(filename=paste(outprefix, ".ram.pdf", sep=""), scale=0.45)


