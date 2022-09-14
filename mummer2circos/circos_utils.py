

def get_contigs_coords(record_list):
    from Bio import SeqIO
    import re

    contigs_add = {}


    contig_start = 0
    contig_end = 0
    for record in record_list:
        name = re.sub("\|", "", record.name)
        contig_start = contig_end+1
        contig_end = contig_start+len(record)
        contigs_add[name] = [contig_start, contig_end]
    return contigs_add


def get_karyotype_from_gbk_or_fasta_single_ref(reference_records,
                                     filter_list=[],
                                     out="circos_contigs.txt",
                                     format='genbank',
                                     concatenated_coordinates=True):
    '''

    filter list = list of contig names to removes: default: empty list
    format = either genbank of fasta
    concatenated coordinates: either cumulate the length of each contigs or use distinct coordinates for each contig

    RETURN: names of the first and last records (usefull to increases intervall between contigs in circos)

    '''


    from Bio import SeqIO
    import re




    contig2start_psoition = {}

    with open(out, 'w') as f:
        # chr - Rhab Rhab 0 1879212 spectral-5-div-4
        i = 0
        contig_start = 0
        contig_end = 0
        for n, record in enumerate(reference_records):
            if n==0:
                first_record = record.name
            current_name = record.name
            if '|' in current_name:
                raise('invalid contig name:\n%s\n' % current_name)
            if current_name in filter_list:
                continue
            # cumulated length if not link plot and not filter_ref (TODO, put it as an argument?)
            if concatenated_coordinates:
                contig_start = contig_end+1
            contig_end = contig_start+len(record)

            # keep in memory
            contig2start_psoition[current_name] = contig_end+1

            # for color scale
            if i == 4:
                i =0
            i+=1

            f.write("chr - %s %s %s %s greys-3-seq-%s\n" % (current_name, current_name, contig_start, contig_end, i)) # i,



    return (first_record, current_name)



def records2locus_tag2start_stop(records):


    '''
    RETURN A DICTIONNARY OF DICTIONNARIES: one/record
    Contains coordinates of each locus
    '''

    record2locus_tag2start_stop = {}

    for record in records:
        record2locus_tag2start_stop[record.name] = {}
        for feature in record.features:
            if feature.type == 'CDS':
                record2locus_tag2start_stop[record.name][feature.qualifiers['locus_tag'][0]] = [int(feature.location.start), int(feature.location.end)]
    return records2locus_tag2start_stop


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def parse_smatools_depth(samtools_depth_file):
    import pandas
    with open(samtools_depth_file, 'r') as f:
        table = pandas.read_csv(f, sep='\t', header=None, index_col=0)
    return table

def get_median_contig_coverage(contig_name_list, samtools_depth_file):
    import numpy
    table = parse_smatools_depth(samtools_depth_file)
    contig2median_cov = {}
    for index in contig_name_list:
        subset_table = table.loc[index]
        contig2median_cov[index] = numpy.median(subset_table.iloc[:,1])
    contig2median_cov['assembly'] = numpy.median(table.iloc[:,1])

    return contig2median_cov


def contig2gc(records):
    from Bio.SeqUtils import GC
    contig2gc_dico = {}
    for record in records:
        contig2gc_dico[record.name] = GC(record.seq)
    return contig2gc_dico


def macsyfinder_table2circos(macsyfinder_table, gbk_record, contigs_add, outfile_name="circos_macsyfinder.txt"):
    from Bio import SeqRecord

    # keep:
    # locus_tag
    # gene
    # system
    # status (mandatory or not)
    hit_data = []



    locus2location = {}
    if type(gbk_record) == list:

        for record in gbk_record:

            for feature in record.features:
                if feature.type == 'CDS':

                    locus2location[feature.qualifiers['locus_tag'][0]] = [feature.location, record.name]
    elif isinstance(gbk_record, SeqRecord):
        for feature in gbk_record.features:
            if feature.type == 'CDS':
                locus2location[feature.qualifiers['locus_tag'][0]] = [feature.location, gbk_record.name]
    else:
        raise IOError('wrong input format')

    all_systems = []
    all_labels=[]


    with open(macsyfinder_table, 'r') as f:
        for i, line in enumerate(f):
            if i==0:
                continue
            else:

                data = line.rstrip().split()

                if len(data) == 16:
                    shift=1
                else:
                    shift=0
                gene = data[5-shift]
                locus_tag= data[0]
                system=data[7]
                if system not in all_systems:
                    all_systems.append(system)
                status=data[10]
                record_name = locus2location[locus_tag][1]
                start = str(locus2location[locus_tag][0].start+contigs_add[record_name][0])
                end = str(locus2location[locus_tag][0].end+contigs_add[record_name][0])
                if status == 'forbidden':
                    continue
                elif status == 'mandatory':
                    label = '%s_M' % (gene)
                else:
                    label = '%s' % (gene)

                all_labels.append([record_name, start, end, label, system])

    c=1
    system2color={}
    for system in all_systems:
        if c==9:
            c=1
        color = 'set1-8-qual-%s' % c
        system2color[system] = color
        c+=1
    o = open(outfile_name, 'w')
    for label in all_labels:
        o.write("%s\t%s\t%s\t%s\tcolor=%s\n" % (label[0], label[1], label[2], label[3], system2color[label[4]]))
    o.close()


def samtools_depth2circos_data(samtools_depth_file, contigs_add, i, window=1000):
    import numpy
    contig2coverage = {}
    all_positions_coverage = []
    with open(samtools_depth_file, 'r') as f:
        for line in f:
            data = line.rstrip().split('\t')
            if data[0] not in contig2coverage:
                contig2coverage[data[0]] = []
                contig2coverage[data[0]].append(int(data[2]))
            else:
                contig2coverage[data[0]].append(int(data[2]))
            all_positions_coverage.append(int(data[2]))
    all_contigs_median = float(numpy.median(all_positions_coverage))
    with open('circos_samtools_depth_%s.txt' % i, 'w') as g:
        for contig in contig2coverage:
            # split list by chunks of 1000
            mychunks = chunks(contig2coverage[contig], window)
            for i, cov_list in enumerate(mychunks):
                median_depth = numpy.median(cov_list)
                if median_depth > (4*all_contigs_median):
                    median_depth = 4*all_contigs_median
                if contigs_add:
                    g.write("%s\t%s\t%s\t%s\n" % (contig, (i*window)+contigs_add[contig][0],
                                                  ((i*window)+(window-1))+contigs_add[contig][0],
                                                  median_depth))
                else:
                    g.write("%s\t%s\t%s\t%s\n" % (contig, (i*window),
                                                  ((i*window)+(window-1)),
                                                  median_depth))
    return all_contigs_median

