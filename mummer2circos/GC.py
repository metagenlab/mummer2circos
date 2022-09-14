#!/usr/bin/env python

# Perform various GC calculations: GC skew, gc variations from the average, plots of length vs gc, circos input files,...
# TODO separate circos stuff from the rest
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

from Bio.SeqUtils import GC, GC_skew

def circos_gc_var(record, windows=1000, shift=0):
    '''
    :param record:
    :return: circos string with difference as compared to the average GC
    ex: average = 32
        GC(seq[3000:4000]) = 44
        diff = 44 - 32 = 12%
    '''
    circos_string = ''
    from Bio.SeqFeature import FeatureLocation
    average_gc = GC(record.seq)
    gap_locations = []
    for feature in record.features:
        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        gap_locations.append(FeatureLocation(0, len(record.seq)))
    else:
        gap_locations.append(FeatureLocation(gap_locations[-1].end + 1, len(record.seq)))
    if len(gap_locations) > 1:
        #gap_locations.append(FeatureLocation(gap_locations[-1].end + 1, len(record.seq)))


        for i in range(0, len(gap_locations)):
            if i == 0:
                seq = record.seq[0:gap_locations[i].start]
                chr_start = 0
            else:
                seq = record.seq[gap_locations[i-1].end:gap_locations[i].start]
                chr_start = gap_locations[i-1].end
            contig_name = record.name + "_%s" % (i +1)

            if len(seq) < windows:
                windows_range = len(seq)
            else:
                windows_range = windows
            if len(seq) == 0:
                continue
            for i in range(0, len(seq), windows_range):
                start = i
                stop = i + windows
                if 'n' in record.seq[start:stop]:
                    continue
                if 'N' in record.seq[start:stop]:
                    continue
                gc = GC(record.seq[start:stop]) - average_gc
                if stop > len(seq):
                    stop = len(seq)
                section_start = chr_start + start
                section_end = chr_start + stop
                circos_string += "%s %s %s %s\n" % (contig_name, section_start, section_end, gc)
    else:
        seq = record.seq
        contig_name = record.id #.split('.')[0]
        for i in range(0, len(seq), windows):
            start = i
            stop = i + windows
            gc = GC(record.seq[start:stop]) - average_gc
            if stop > len(seq):
                stop = len(seq)
                if stop - start < 500:
                    break
            circos_string += "%s %s %s %s\n" % (contig_name, start+shift, stop+shift, gc)
    return circos_string



def circos_gc_content(record, windows=1000, shift=0):
    '''
    :param record:
    :return: circos string with difference as compared to the average GC
    ex: average = 32
        GC(seq[3000:4000]) = 44
        diff = 44 - 32 = 12%

    UPDATE 12.06.2017: calculs based on complete (concatenated) sequence, then converted to draft contigs coords
    '''
    circos_string = ''
    from Bio.SeqFeature import FeatureLocation
    average_gc = GC(record.seq)
    gap_locations = []
    for feature in record.features:
        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        gap_locations.append(FeatureLocation(0, len(record.seq)))
    else:
        gap_locations.append(FeatureLocation(gap_locations[-1].end + 1, len(record.seq)))
    if len(gap_locations) > 1:
        #gap_locations.append(FeatureLocation(gap_locations[-1].end + 1, len(record.seq)))




        for i in range(0, len(gap_locations)):
            if i == 0:
                seq = record.seq[0:gap_locations[i].start]
                chr_start = 0
            else:
                seq = record.seq[gap_locations[i-1].end:gap_locations[i].start]
                chr_start = gap_locations[i-1].end
            contig_name = record.name + "_%s" % (i +1)
            if len(seq) <= windows:
                window_range = len(seq)
            else:
                window_range = windows
            for i in range(0, len(seq), window_range):
                start = i
                stop = i + windows
                gc = GC(record.seq[start:stop])
                if stop > len(seq):

                    stop = len(seq)
                section_start = chr_start + start
                section_end = chr_start + stop
                circos_string += "%s %s %s %s\n" % (contig_name, section_start, section_end, gc)

    else:
        seq = record.seq
        contig_name = record.id #.split('.')[0]
        for i in range(0, len(seq), windows):
            start = i
            stop = i + windows
            gc = GC(record.seq[start:stop])
            if stop > len(seq):
                stop = len(seq)
                if stop - start < 500:
                    break
            circos_string += "%s %s %s %s\n" % (contig_name, start+shift, stop+shift, gc)
    return circos_string



def circos_gc_skew(record, windows=1000, shift=0):
    '''
    :param record:
    :return: circos string with difference as compared to the average GC
    ex: average = 32
        GC(seq[3000:4000]) = 44
        diff = 44 - 32 = 12%

    NEW 12.06.2017: calculate GC based on whole sequence

    '''
    from Bio.SeqFeature import FeatureLocation
    circos_string = ''

    gap_locations = []
    for feature in record.features:

        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        gap_locations.append(FeatureLocation(0, len(record.seq)))

    else:

        gap_locations.append(FeatureLocation(len(record.seq), len(record.seq)))
    if len(gap_locations) > 1:
        from mummer2circos import circos_convert_contigs_coords

        contig_coords =  []
        start = 0
        for i, coord in enumerate(gap_locations):
            contig_name = record.name + "_%s" % (i + 1)
            contig_coords.append([contig_name, start, int(coord.start)])
            start = coord.end + 1

        values = GC_skew(record.seq, windows)
        data_list = []
        for n_value, value in enumerate(values):

            start = (windows * n_value) + 1
            end = (start + windows)
            #for i in range(0, len(gap_locations)):
            data_list.append([record.name, start, end, value])

        renamed_data = circos_convert_contigs_coords.rename_karyotype(contig_coords, data_list)
        for row in renamed_data:
            contig_name = row[0]
            start = row[1]
            end = row[2]
            try:
                value = row[3]
            except:
                continue
            circos_string += "%s %s %s %s\n" % (contig_name, start, end, value)

    else:
        try:
            values = GC_skew(record.seq, windows)
        except:
            values = GC_skew(record.seq, 2001)
        # skip last value
        for i in range(0, len(values)-1):

            start = i *windows
            stop = start + windows

            circos_string += "%s %s %s %s\n" % (record.id, start+shift, stop+shift, values[i])
    return circos_string


def circos_cumul_gc_skew(record, windows=1000, shift=0, initial=0):
    '''
    :param record:
    :return: circos string with difference as compared to the average GC
    ex: average = 32
        GC(seq[3000:4000]) = 44
        diff = 44 - 32 = 12%

    '''
    from Bio.SeqFeature import FeatureLocation
    circos_string = ''

    gap_locations = []
    for feature in record.features:

        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        gap_locations.append(FeatureLocation(0, len(record.seq)))

    else:
        gap_locations.append(FeatureLocation(len(record.seq), len(record.seq)))

    if len(gap_locations) > 1:
        for i in range(0, len(gap_locations)):
            if i == 0:
                seq = record.seq[0:gap_locations[i].start]
                chr_start = 0
            else:
                seq = record.seq[gap_locations[i-1].end:gap_locations[i].start]
                chr_start = gap_locations[i-1].end
            try:
                values = GC_skew(seq, windows)
            except:
                pass
            contig_name = record.name + "_%s" % (i + 1)

            for i in range(0, len(values)):
                start = i *windows
                stop = start + windows
                section_start = chr_start + start
                section_end = chr_start + stop
                circos_string += "%s %s %s %s\n" % (contig_name, section_start+shift, section_end+shift, values[i])
    else:
        try:
            values = GC_skew(record.seq, windows)
        except:
            try:
                values = GC_skew(record.seq, 2001)
            except:
                print(record.seq)
                import sys
                sys.exit()
        # skip last value
        acc = initial
        for i in range(0, len(values)-1):
            gc = values[i]
            acc += gc
            start = i * windows
            stop = start + windows

            circos_string += "%s %s %s %s\n" % (record.id, start+shift, stop+shift, acc)
    return circos_string, acc









