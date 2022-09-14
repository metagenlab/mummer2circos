#!/usr/bin/env python
# python 2.7.5 requires biopython

########### promer2circos ############

def get_contig(location, contig_coordlist, contig_pos):
    for i, one_contig in enumerate(contig_coordlist):
        if location >= one_contig[1] and location <= one_contig[2]:
            return one_contig, contig_coordlist[i + 1:len(contig_coordlist)], location
        elif location > contig_coordlist[i - 1][2] and location <= one_contig[1]:
            if contig_pos == 'start':
                return one_contig, contig_coordlist[i + 1:len(contig_coordlist)], one_contig[1]
            else:
                return contig_coordlist[i - 1], contig_coordlist[i:len(contig_coordlist)], contig_coordlist[i - 1][2]

        else:
            # partial match, probably overlap between end contig and gap region
            continue
    return False, False, False


def rename_karyotype(contig_coords, data_list):
    '''

    :param contig_coords: chr - NZ_JSAM00000000_1 NZ_JSAM00000000_1 0 104228 spectral-5-div-4 ==> keep 3, 4 et 5
    :param data_list: list of contig coords: [[contig_X, start, end],[...]]
    :return:
    '''

    import copy

    renamed_data = []
    for i, data in enumerate(data_list):
        start = int(data[1])
        end = int(data[2])
        contig_start, following_contigs1, position1 = get_contig(start, contig_coords, contig_pos="start")

        contig_end, following_contigs2, position2 = get_contig(end, contig_coords, contig_pos="end")

        if (contig_start is False) or (contig_end is False):
            if start > contig_coords[-1][1]:
                contig_start, following_contigs1, position1 = contig_coords[-1], [], start
                contig_end, following_contigs2, position2 = contig_coords[-1], [], contig_coords[-1][2]
        if contig_end is False:
            continue
        data[1] = position1
        data[2] = position2

        if contig_start[0] == contig_end[0]:
            data[0] = contig_start[0]
            renamed_data.append(data)
        else:
            # span across 2 contigs: make 2 coordinates (until the end of the first contig and from the begining of the second)
            data_1 = copy.copy(data)
            data_1[0] = contig_start[0]
            data_1[2] = contig_start[2]
            renamed_data.append(data_1)
            # enumerate following contigs until we match the final one
            for contig2 in following_contigs1:
                # final contig of the range, add it and break the inner loop
                if contig2[0] == contig_end[0]:
                    data_2 = copy.copy(data)
                    data_2[0] = contig_end[0]
                    # start from the first position of the second contig
                    data_2[1] = contig_end[1]
                    renamed_data.append(data_2)
                    break
                else:
                    # entire contig comprised within the range
                    # add it entiely to the new list
                    renamed_data.append(contig2)


    return renamed_data


def read_circos_file(circos_file):
    data_list = []
    with open(circos_file) as f:
        for row in f:
            data = row.rstrip().split(' ')
            if len(data) < 3:
                data = row.rstrip().split('\t')
            data_list.append(data)

    return data_list

