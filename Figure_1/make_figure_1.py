#!/usr/bin/env python3

# Makes Figure 1 of the paper

import os
import math
import fastaq
import iva


iva_contig_name = 'contig.00001'
iva_fa = 'contigs.iva.fa'
price_fa = 'contigs.price.fa'
trinity_fa = 'contigs.trinity.fa'
vicuna_fa = 'contigs.vicuna.fa'
iva_read_depth = 'contigs.iva.fa.read_depth.gz'
outprefix = 'figure_1'


class Coords:
    def __init__(self, a, b):
        self.start = a
        self.end = b

    def __lt__(self, other):
        return self.start < other.start or (self.start == other.start and self.end > other.end)

    def __str__(self):
        return str(self.start) + ',' + str(self.end)


def nucmer_to_coords(refname, filename):
    f = fastaq.utils.open_file_read(filename)
    lines = f.readlines()
    fastaq.utils.close(f)
    i = 0
    while not lines[i].startswith('[S1]'):
        i += 1
    lines = [x.rstrip() for x in lines[i+1:]]
    coords = []

    for line in lines:
        a = line.split('\t')
        if a[11] == refname:
            s = min(int(a[0]), int(a[1]))
            e = max(int(a[0]), int(a[1]))
            coords.append(Coords(s, e))
    coords.sort()
    return coords


def coords_to_plot_positions(coords):
    pos = []

    for c in coords:
        if len(pos) == 0:
            pos = [[c]]
            continue

        possible = [i for i in range(len(pos)) if c.start > pos[i][-1].end + 1]
        if len(possible) == 0:
            pos.append([c])
        else:
            first_i = possible[0]
            for  i in possible:
                if pos[i][-1].end < pos[first_i][-1].end:
                    first_i = i

            pos[first_i].append(c)

    return pos


def load_cov_plot(filename, refname):
    lines = fastaq.utils.syscall_get_stdout('tabix ' + filename + ' ' + refname)
    read_depth = []
    snps = []
    for line in lines:
        a = line.split()
        read_depth.append(int(a[2]))
        snps.append(int(a[3]))

    return read_depth, snps


def list_to_plot_coords(l, xmin, xmax, ymin, ymax, y_axis_max, window=1):
    max_y = max(l)
    min_y = min(l)
    coords = []
    for i in range(0, len(l), window):
        y = sum(l[i:i+window]) / window
        y = 1 - (y / y_axis_max) # scaled between 0 and 1
        y = ymin + (y * (ymax - ymin)) # scaled to plot area
        x = i / len(l) # scaled between 0 and 1
        x = xmin + (x * (xmax - xmin))
        coords.append((x, y))

    return coords


# coords  = list of tuples [(x1, y1), (x2, y2) ...]
def svg_polygon(coords, fill_colour, border_colour, border_width = 1, opacity=-1):
    return_string = '<polygon points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
            + 'fill="' + fill_colour + '" '

    if opacity != -1:
        return_string += 'fill-opacity="' + str(opacity) + '" '

    return_string += 'stroke="' + border_colour + '" ' \
            + 'stroke-width="' + str(border_width) + '" ' \
            + '/>'

    return return_string


def svg_polyline(coords, colour, thickness = 1):
    return '<polyline points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
            + 'stroke="' + colour + '" ' \
            + 'stroke-width="' + str(thickness) + '" ' \
            + 'fill="none"' \
            + '/>'


def svg_text(x, y, text, vertical=False, size=12):
    s = '<text x="' + str(x) + '" y="'  + str(y) + '" font-size="' + str(size) + '"'
    if vertical:
        s += ' transform="rotate(-90,' + str(x) + ',' + str(y) + ')"'
    s += '>' + text + '</text>'
    return s


coords = {}
vicuna_coords_file = outprefix + '.vicuna.coords'
price_coords_file = outprefix + '.price.coords'
trinity_coords_file = outprefix + '.trinity.coords'
iva.mummer.run_nucmer(vicuna_fa, iva_fa, vicuna_coords_file, min_id=80, min_length=100, breaklen=500)
iva.mummer.run_nucmer(price_fa, iva_fa, price_coords_file, min_id=80, min_length=100, breaklen=500)
iva.mummer.run_nucmer(trinity_fa, iva_fa, trinity_coords_file, min_id=80, min_length=100, breaklen=500)
coords['VICUNA'] = nucmer_to_coords(iva_contig_name, vicuna_coords_file)
coords['PRICE'] = nucmer_to_coords(iva_contig_name, price_coords_file)
coords['Trinity'] = nucmer_to_coords(iva_contig_name, trinity_coords_file)


iva_lengths = {}
fastaq.tasks.lengths_from_fai(iva_fa + '.fai', iva_lengths)
contig_height = 2
contig_space = 2
assembly_space = 4
contig_rows = sum([len(x) for x in coords.values()])
total_contig_height = contig_rows * (contig_height + contig_space) + len(coords.values()) * assembly_space
plot_height = 50
plot_space = 20
gap_iva_to_axis = 3
height = total_contig_height + plot_height * 2 + plot_space * 2 * gap_iva_to_axis
axis_thickness = 0.5


f = fastaq.utils.open_file_write(outprefix + '.svg')
total_width = 520
plotting_width = 450
plot_x_start = total_width - plotting_width
print (r'''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg width="''' + str(total_width) + '" height="' + str(height) + '">', file=f)

colours = {
    'PRICE': 'white',
    'Trinity': 'black',
    'VICUNA': 'white',
}

h = height
for assembler in ['VICUNA', 'Trinity', 'PRICE']:
    positions = coords_to_plot_positions(coords[assembler])
    for l in positions:
        top = h - contig_height
        bottom = h
        for p in l:
            left = plot_x_start + plotting_width * p.start / iva_lengths[iva_contig_name]
            right = plot_x_start + plotting_width * p.end / iva_lengths[iva_contig_name]
            corners = [(left, bottom), (left, top), (right, top), (right, bottom)]
            print(svg_polygon(corners, colours[assembler], 'black', border_width=0.25), file=f)
        h -= (contig_height + contig_space)
    print(svg_text(10, h + assembly_space / 2, assembler, size=12), file=f)
    h -= assembly_space


top = h - contig_height
bottom = h
print(svg_polygon([(plot_x_start, bottom), (plot_x_start, top), (total_width, top), (total_width, bottom)], 'black', 'black', border_width=0.25), file=f)
print(svg_text(10, h, "IVA", size=12), file=f)

h -= plot_space

print(svg_polyline([(plot_x_start,h), (total_width, h)], 'black', thickness=axis_thickness), file=f) # x axis
# x axis ticks and labels
iva_length = iva_lengths[iva_contig_name]
print(svg_text(plot_x_start-2, h+9, '0', size=10), file=f)
print(svg_polyline([(plot_x_start,h), (plot_x_start,h+3)], 'black', thickness=axis_thickness), file=f)
for i in range(1,9,2):
    tick_x = 1000 * (i + 1)
    x_pos = plot_x_start + (tick_x / iva_length) * plotting_width
    print(svg_text(x_pos-11, h+9, str(tick_x), size=10), file=f) 
    print(svg_polyline([(x_pos,h), (x_pos,h+3)], 'black', thickness=axis_thickness), file=f)

h -= gap_iva_to_axis


read_depth, snps = load_cov_plot(iva_read_depth, iva_contig_name)
max_read_depth = round(max(read_depth), -3)
read_depth_coords = list_to_plot_coords(read_depth, plot_x_start, total_width, h - plot_height, h, max_read_depth, window=1)
print(svg_polyline(read_depth_coords, 'black', thickness=0.25), file=f)
print(svg_polyline([(plot_x_start,h), (total_width, h)], 'gray', thickness=0.1), file=f)
print(svg_polyline([(plot_x_start - 5, h), (plot_x_start - 5, h - plot_height)], 'black', thickness=axis_thickness), file=f) # y axis
print(svg_polyline([(plot_x_start - 8, h), (plot_x_start - 5, h)], 'black', thickness=axis_thickness), file=f) # tick
print(svg_polyline([(plot_x_start - 8, h - plot_height), (plot_x_start - 5, h - plot_height)], 'black', thickness=axis_thickness), file=f) # tick


print(svg_text(plot_x_start - 15, h, "0", size=10), file=f)  # tick label
print(svg_text(plot_x_start - 29, h - plot_height, str(max_read_depth), size=10), file=f)  # tick label
print(svg_text(10, h, "Read depth", size=12, vertical=True), file=f)

h -= plot_space + plot_height
snps_per_read = [snps[i] / read_depth[i] if read_depth[i] > 0 else 0 for i in range(len(snps))]
snps_per_read_coords =  list_to_plot_coords(snps_per_read, plot_x_start, total_width, h - plot_height, h, 1, window=1)
print(svg_polyline(snps_per_read_coords, 'black', thickness=0.25), file=f)
print(svg_polyline([(plot_x_start,h), (total_width, h)], 'gray', thickness=0.1), file=f)
print(svg_polyline([(plot_x_start-5,h), (plot_x_start-5, h - plot_height)], 'black', thickness=axis_thickness), file=f) # y axis
print(svg_polyline([(plot_x_start - 8, h), (plot_x_start - 5, h)], 'black', thickness=axis_thickness), file=f) # tick
print(svg_polyline([(plot_x_start - 8, h - plot_height), (plot_x_start - 5, h - plot_height)], 'black', thickness=axis_thickness), file=f) # tick
print(svg_text(plot_x_start-15, h, "0", size=10), file=f)  # tick label
print(svg_text(plot_x_start - 15, h - plot_height, "1", size=10), file=f)  # tick label
print(svg_text(10, h, "% read differences", size=12, vertical=True), file=f)

print('</svg>', file=f)
fastaq.utils.close(f)


os.unlink(outprefix + '.price.coords')
os.unlink(outprefix + '.trinity.coords')
os.unlink(outprefix + '.vicuna.coords')
