#!/usr/bin/env python

########### promer2circos ############

# A script to generate circos plots from one reference genbank and fasta files
# trestan.pillonel@gmail.com

######################################

import sys;
import re;
from mummer2circos import circos_utils
import random
import string
import logging
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s',datefmt='%H:%M:%S')

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def purge(dir, pattern):
    import os
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


class CircosException(Exception):
    pass


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


class Fasta2circos():
    def __init__(self, fasta1, fasta2,
                 filter_ref=True,
                 filter_query=True,
                 heatmap=False,
                 samtools_depth=False,
                 gaps=False,
                 blast=False,
                 gc=True,
                 highlight_list=[],
                 algo="nucmer",
                 min_gap_size=1000,
                 blastn=False,
                 gbk2orf=False,
                 window_size=1000,
                 secretion_systems=False,
                 condensed_tracks=False,
                 label_file=False,
                 locus_taxonomy=False,
                 blast_identity_cutoff=80,
                 force_data_dir=False):

        from mummer2circos import nucmer_utility
        import os
        from Bio import SeqIO
        from mummer2circos import gbk2circos
        import shutil

        self.contigs_add = {}
        self.working_dir = os.getcwd()
        self.last_track = 0.99
        self.window_size = window_size
        self.locus2taxo_color = None
        self.gbk2orf = gbk2orf
        self.blast = blast
        self.heatmap = heatmap
        self.samtools_depth = samtools_depth
        self.circos_data_dir = "circos_data"
        if os.path.exists(self.circos_data_dir):
            if not force_data_dir:
                raise IOError("circos data folder already present in the working directory, remove it or use -f to force rewrite")
            else:
                shutil.rmtree(self.circos_data_dir)
        if os.path.exists("circos_tmp"):
                shutil.rmtree("circos_tmp")
        os.mkdir(self.circos_data_dir)
        os.mkdir("circos_tmp")
        
        logging.info(f"Genomes will be aligned with:\t{algo}")
        
        if locus_taxonomy:
            self.parse_taxonomy_file(locus_taxonomy)
        if fasta2 is not None:
            if algo == "nucmer":
                logging.info(f"Reference genome:\t{fasta1}")
                nucmer_utility.execute_promer(fasta1, fasta2, algo="nucmer", out_dir="circos_tmp")
            elif algo == "megablast":
                self.execute_megablast(fasta1, fasta2)
            elif algo == "promer":
                nucmer_utility.execute_promer(fasta1, fasta2, algo="promer", out_dir="circos_tmp")
        else:
             fasta2=[]
        if not heatmap:
                outname = os.path.basename(fasta2[0]).split('.')[0]
                hit_list, query_list = self.get_link("circos_tmp/%s.coords" % outname, 
                                                     link_file=f"{self.circos_data_dir}/circos.link",
                                                     algo=algo)
        else:

            # coords for cumulated coorginates
            records = [i for i in SeqIO.parse(open(fasta1), "fasta")]
            self.contigs_add = circos_utils.get_contigs_coords(records)

            self.circos_reference = gbk2circos.Circos_config(f"{self.circos_data_dir}/circos_contigs.txt",
                                                             show_ideogram_labels="no",
                                                             radius=0.7,
                                                             show_tick_labels="yes",
                                                             show_ticks="yes")

            if gbk2orf:
                minus_file, plus_file = self.gbk2circos_data(gbk2orf, minus_file=f"{self.circos_data_dir}/circos_orf_minus.txt", plus_file=f"{self.circos_data_dir}/circos_orf_plus.txt")
                
                self.circos_reference.add_highlight(os.path.basename(minus_file),
                                                    'grey_a1',
                                                    r1="%sr" % (self.last_track - 0.01),
                                                    r0="%sr" % (self.last_track - 0.03))
                self.last_track -= 0.03

                self.circos_reference.add_highlight(os.path.basename(plus_file),
                                                    'grey_a1',
                                                    r1="%sr" % (self.last_track),
                                                    r0="%sr" % (self.last_track - 0.02))
                self.last_track -= 0.05

            # genome_list = [i.split('.')[0] + '.heat' for i in fasta2]

            # self.add_multiple_genome_tracks(genome_list, highlight_list)

            updated_list = []
            gap_hilight_list = []
            all_hit_list = []
            all_query_list = []
            for i, one_fasta in enumerate(fasta2):
                logging.info(f"Query genome {i+1}/{len(fasta2)}: {one_fasta}")
                out_prefix = os.path.basename(one_fasta).split('.')[0]

                if i % 2 == 0:
                    col = 200
                else:
                    col = 250
                # try:
                if algo == "nucmer" or algo == 'promer':
                    # hit_list, query_list, contig2start_stop_list = self.nucmer_coords2heatmap("%s.coords" % one_fasta.split('.')[0], col=col, algo=algo)

                    hit_list, query_list = nucmer_utility.coord_file2circos_heat_file("circos_tmp/%s.coords" % out_prefix,
                                                                                      self.contigs_add, 
                                                                                      algo=algo,
                                                                                      out_dir=self.circos_data_dir)
                    if gaps:
                        contig2start_stop_list, coords_file = nucmer_utility.delta_file2start_stop_list("circos_tmp/%s.delta" % out_prefix,
                                                                                                        self.contigs_add, 
                                                                                                        algo=algo,
                                                                                                        minimum_identity=5,
                                                                                                        out_dir="circos_tmp")

                elif algo == "megablast":
                    hit_list, query_list, contig2start_stop_list = self.megablast2heatmap("circos_tmp/blast_result_%s.tab" % out_prefix, col=col)
                else:
                    raise IOError('unknown algo!')
                all_hit_list += hit_list
                all_query_list += query_list
                updated_list.append(os.path.basename(one_fasta))
                # except:
                #    continue
                if gaps:

                    nucmer_utility.get_circos_gap_file_from_gap_start_stop(contig2start_stop_list,
                                                                           self.contigs_add,
                                                                           out_highlight=f"{self.circos_data_dir}/circos_gaps_highlight_%s.txt" % i,
                                                                           out_labels=f"{self.circos_data_dir}/circos_gaps_labels_%s.txt" % i,
                                                                           min_gap_size=min_gap_size,
                                                                           gap_merge_distance=0)
                    # get_gaps_from_start_stop_lists(contig2start_stop_list, out_highlight="circos_gaps_highlight_%s.txt" % i, out_labels="circos_gaps_labels_%s.txt" % i, min_gap_size=min_gap_size)
                    if i % 2 == 0:
                        color = "orrd-9-seq-9"
                    else:
                        color = "blues-9-seq-9"
                
                    gap_hilight_list.append([color, f"{self.circos_data_dir}/circos_gaps_highlight_%s.txt" % i])

                    supp = '''show_links     = yes
                            link_dims      = 4p,4p,4p,4p,4p
                            link_thickness = 2p
                            link_color     = red

                            label_size   = 24p
                            label_font   = condensed

                            padding  = 0p
                            rpadding = 0p
                            '''
                    # self.circos_reference.add_plot("circos_gaps_labels_%s.txt" % i ,type="text", r0="1r", r1="1.3r",color="black",rules=supp)
            self.genome_list = [f'{i.split(".")[0]}.heat' for i in updated_list]
            self.add_multiple_genome_tracks(self.genome_list, highlight_list, condensed=condensed_tracks)
            if gaps:
                for i in gap_hilight_list:
                    self.circos_reference.add_highlight(i[1],
                                                        i[0],
                                                        r1="%sr" % (self.last_track - 0.01),
                                                        r0="%sr" % (self.last_track - 0.03))
                    self.last_track -= 0.03

            if blast:
                self.blast2circos_file(blast, fasta1, blastn=blastn, identity_cutoff=blast_identity_cutoff, out_folder=self.circos_data_dir)

                # self.circos_reference.add_highlight("circos_blast.txt",
                #                               'red',
                #                               r1="%sr" % (self.last_track - 0.01),
                #                               r0="%sr" % (self.last_track - 0.03))
                # self.last_track -= 0.03
                # snuggle_refine                 = yes
                supp = '''
                        label_snuggle             = yes

                        max_snuggle_distance            = 20r
                        show_links     = yes
                        link_dims      = 10p,88p,30p,4p,4p
                        link_thickness = 2p
                        link_color     = red

                        label_size   = 24p
                        label_font   = condensed

                        padding  = 0p
                        rpadding = 0p
                        '''
                self.circos_reference.add_plot(f"{self.circos_data_dir}/circos_blast_labels.txt", 
                                               type="text", 
                                               r0="1r", 
                                               r1="2r", 
                                               color="black",
                                               rules=supp)

            if label_file:
                if blast:
                    dist="288"
                else:
                    dist="120"

                supp = '''
                        label_snuggle             = yes

                        max_snuggle_distance            = 20r
                        show_links     = yes
                        link_dims      = 10p,%sp,30p,4p,4p
                        link_thickness = 2p
                        link_color     = blue

                        label_size   = 24p
                        label_font   = condensed

                        padding  = 0p
                        rpadding = 0p
                        ''' % dist

                o = open('circos.labels.tab', 'w')
                with open(label_file, 'r') as f:
                    for line in f:
                        data = line.rstrip().split('\t')
                        start = str(int(data[1]) + self.contigs_add[data[0]][0])
                        end = str(int(data[2]) + self.contigs_add[data[0]][0])
                        data[1] = start
                        data[2] = end
                        o.write('\t'.join(data) + '\n')
                o.close()

                self.circos_reference.add_plot(f'{self.circos_data_dir}/circos.labels.tab', type="text", r0="1r", r1="2r", color="black",
                                               rules=supp)

            if gbk2orf and secretion_systems:
                self.records = [i for i in SeqIO.parse(gbk2orf, 'genbank')]
                circos_utils.macsyfinder_table2circos(secretion_systems, self.records, self.contigs_add,
                                                      'circos_secretion_systems.txt')

                supp = '''
                        label_snuggle             = yes

                        max_snuggle_distance            = 20r
                        show_links     = yes
                        link_dims      = 10p,128p,30p,4p,4p
                        link_thickness = 2p
                        link_color     = red
                        snuggle_link_overlap_test = yes
                        snuggle_link_overlap_tolerance = 10p
                        label_size   = 14p
                        label_font   = condensed
                        padding  = 0p
                        rpadding = 0p

                        '''
                self.circos_reference.add_plot('circos_secretion_systems.txt', type="text", r0="1r", r1="1.6r",
                                               color="black", rules=supp)

            if gc:
                from mummer2circos import GC
                from Bio import SeqIO

                fasta_records = list(SeqIO.parse(fasta1, 'fasta'))

                out_var_file = (f'{self.circos_data_dir}/circos_GC_var.txt')
                out_skew_file = (f'{self.circos_data_dir}/circos_GC_skew.txt')
                f = open(out_var_file, 'w')
                g = open(out_skew_file, 'w')

                out_var = ''
                out_skew = ''
                out_cumul_skew = ''
                out_gc_content = ''
                initial = 0
                for n, record in enumerate(fasta_records):
                    contig = re.sub("\|", "", record.id)
                    shift = self.contigs_add[contig][0]
                    # this function handle scaffolds (split sequence when encountering NNNNN regions)
                    out_var += GC.circos_gc_var(record, self.window_size, shift=shift)
                    circos_cumul, initial = GC.circos_cumul_gc_skew(record, self.window_size, shift=shift,
                                                                    initial=initial)
                    out_cumul_skew += circos_cumul
                    out_skew += GC.circos_gc_skew(record, self.window_size, shift=shift)
                    out_gc_content += GC.circos_gc_content(record, self.window_size, shift=shift)
                f.write(out_var)
                g.write(out_skew)
                f.close()
                g.close()

                rule = """<rule>
                        condition          = var(value) < 0
                        fill_color         = lred
                        color = red
                        </rule>

                        <rule>
                        condition          = var(value) > 0
                        fill_color         = lblue
                        color = blue
                        </rule>
                """

                rule2 = """<rule>
                        condition          = var(value) < 0
                        fill_color         = lgreen
                        color = green
                        </rule>

                        <rule>
                        condition          = var(value) > 0
                        fill_color         = lblue
                        color = blue
                        </rule>
                """
                rule3 = """<rule>
                        condition          = var(value) < 0
                        fill_color         = lgreen
                        color = green
                        </rule>

                        <rule>
                        condition          = var(value) > 0
                        fill_color         = lblue
                        color = blue
                        </rule>
                """

                self.last_track = self.last_track - 0.1

                conditions = self.circos_reference.template_rules % (rule)
                self.circos_reference.add_plot(os.path.basename(out_skew_file), 
                                               fill_color="green",
                                               r0="%sr" % (self.last_track - 0.01), r1="%sr" % (self.last_track - 0.09),
                                               type="line", rules=conditions)

                conditions = self.circos_reference.template_rules % (rule2)
                self.circos_reference.add_plot(os.path.basename(out_var_file), 
                                               fill_color="green",
                                               r0="%sr" % (self.last_track - 0.11), r1="%sr" % (self.last_track - 0.19),
                                               type="line", rules=conditions)
                self.last_track = self.last_track - 0.11

                '''
                out_cumul = open('circos_cumul_GC_skew.txt', 'w')
                out_cumul.write(out_cumul_skew)
                out_cumul.close()
                conditions = self.circos_reference.template_rules % (rule3)
                self.circos_reference.add_plot('circos_cumul_GC_skew.txt', fill_color="white", r0="%sr" % (self.last_track -0.11), r1= "%sr" % (self.last_track -0.19), type="line", rules=conditions)
                self.last_track = self.last_track -0.11

                out_cumul = open('circos_GC_content.txt', 'w')
                out_cumul.write(out_gc_content)
                out_cumul.close()
                conditions = self.circos_reference.template_rules % (rule3)
                self.circos_reference.add_plot('circos_GC_content.txt', fill_color="black", r0="%sr" % (self.last_track -0.11), r1= "%sr" % (self.last_track -0.19), type="line", min="0", max="100")
                self.last_track = self.last_track -0.11
                '''

        if heatmap:
            logging.info('Plotting heatmap type plot...')
            self.get_karyotype_from_fasta(fasta1, 
                                          fasta2, 
                                          list(set(all_hit_list)),
                                          list(set(all_query_list)), filter_ref, filter_query,
                                          both_fasta=False,
                                          out=f"{self.circos_data_dir}/circos_contigs.txt")
            # self.config = self.get_circos_config(c1, c2, c3, c4, link=False, heat=True)

            self.config = self.circos_reference.get_file()

            if samtools_depth is not None:
                self.samtools_depth_track_list = []
                for i, depth_file in enumerate(samtools_depth):
                    track_file = f'circos_samtools_depth_%s.txt' % i
                    all_contigs_median = self.samtools_depth2circos_data(depth_file, i, window=self.window_size)
                    self.samtools_depth_track_list.append(track_file)
                    self.add_samtools_depth_track(track_file,
                                                  lower_cutoff=int(all_contigs_median) / 2,
                                                  top_cutoff=int(all_contigs_median) * 2,
                                                  r1=self.last_track - 0.08,
                                                  r0=self.last_track - 0.25)
                    self.last_track -= 0.25

        else:
            logging.info('Plotting link type plot....\n')
            # c1 last_seq_id
            # c2 first_seq_id
            # c3 mid1 last hit id
            # c4 mid2

            last_seq_id, first_seq_id, mid1, mid2 = self.get_karyotype_from_fasta(fasta1, 
                                                                                  fasta2, 
                                                                                  hit_list, 
                                                                                  query_list,
                                                                                  filter_ref, 
                                                                                  filter_query,
                                                                                  both_fasta=True, 
                                                                                  cumul=False,
                                                                                  out=f"{self.circos_data_dir}/circos_contigs.txt")

            self.circos_reference = gbk2circos.Circos_config(f"{self.circos_data_dir}/circos_contigs.txt",
                                                             show_ideogram_labels="yes",
                                                             radius=0.7,
                                                             show_tick_labels="yes",
                                                             show_ticks="yes",
                                                             chr_spacing_list=[[last_seq_id, first_seq_id],
                                                                               [mid1, mid2]],
                                                             ideogram_spacing=0.5,
                                                             color_files='\n<<include colors.my>>')

            self.circos_reference.add_link(f"circos.link")

            if samtools_depth is not None:
                for i, depth_file in enumerate(samtools_depth):
                    all_contigs_median = circos_utils.samtools_depth2circos_data(depth_file, False, i,
                                                                                 window=self.window_size)
                    self.add_samtools_depth_track(f'circos_samtools_depth_%s.txt' % i,
                                                  lower_cutoff=int(all_contigs_median) / 2,
                                                  top_cutoff=int(all_contigs_median) * 2,
                                                  r1=1.45,
                                                  r0=1.25)
                    # self.last_track -= 0.25
            self.config = self.circos_reference.get_file()
            # last_seq_id, first_seq_id, mid1, mid2
            '''
            self.config = self.get_circos_config(last_seq_id,
                                                 first_seq_id,
                                                 mid1,
                                                 mid2, heat=False)
            '''

        self.brewer_conf = """


        # Colors from www.ColorBrewer.org by Cynthia A. Brewer, Geography, Pennsylvania State University.
        # See BREWER for license. See www.colorbrewer.org for details.
        #
        # Color names are formatted as PALETTE-NUMCOLORS-TYPE-COLORCODE
        #
        # where PALETTE is the palette name, NUMCOLORS is the number of colors in the palette, TYPE is the palette type (div, seq, qual) and COLORCODE is the color index within the palette (another versison of the color is defined where COLORCODE is the color's letter, unique to a combination of PALETTE and TYPE)

        #
        # the value after the trailing comment is the palette's critical value. See http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html for details.

        # sequential palettes

        # 5 color palettes

        #violet = 197,27,138
        ortho1 = 199,233,180
        ortho2 = 127,205,187
        ortho3 = 127,205,187
        ref = 254,153,41
        pred = 255,0,55
        pblue = 0, 55, 255
        #highlight = 107, 155, 0
        group_size = 187,255,167
        not_conserved = 254, 29,29
        chlamydiales = 24,116,205
        non_chlamydiales = 0,255,255
        back = 240,240,240
        blue = 1,188,255
        green = 27,255,1
        #sta2 = 255,128,0
        sta1 = 54,144,192
        #euk = 255,131,250

        """

    def parse_taxonomy_file(self, file_path):

        match_color = 'red'
        match_color_2 = 'spectral-5-div-5'
        match_other_color = 'green'
        no_match_color = 'grey_a1'
        match_taxon_2 = ''


        self.locus2taxo_color = {}
        with open(file_path, 'r') as f:
            for n, row in enumerate(f):
                if n == 0:
                    if row[0] == '#':
                        match_taxon = row.rstrip()[1:]
                    else:
                        raise('Taxon rto search for should be indicated as comment on the first line: #Chlamydiae')
                elif n == 1:
                    if row[0] == '#':
                        match_taxon_2 = row.rstrip()[1:]
                else:
                    data = row.rstrip().split('\t')
                    if len(data) == 1:
                        self.locus2taxo_color[data[0]] = no_match_color
                    else:
                        if data[1] == match_taxon:
                            self.locus2taxo_color[data[0]] = match_color
                        elif data[1] == match_taxon_2:
                            self.locus2taxo_color[data[0]] = match_color_2
                        else:
                            self.locus2taxo_color[data[0]] = match_other_color



    def blast2circos_file(self, 
                          blast, 
                          reference, 
                          blastn=False, 
                          identity_cutoff=80,
                          out_folder='.'):

        '''

        tblastn vs contigs by default
        can be switch to blastn

        :param blast:
        :param reference:
        :param blastn:
        :return:
        '''

        from mummer2circos import shell_command
        from mummer2circos import blast_utils
        from Bio.Blast.Applications import NcbitblastnCommandline
        from Bio.Blast.Applications import NcbiblastnCommandline
        import os 
        
        # todo catch IO errors, orther potential errors
        a, b, c = shell_command.shell_command('makeblastdb -in %s -dbtype nucl' % (reference))
        if c != 0:
            raise Exception(f"Problem building blast database:\n{a}\n{b}")
        logging.info(c)
        if not blastn:
            blast_cline = NcbitblastnCommandline(query=blast,
                                                 db=reference,
                                                 evalue=0.00000001,  # 0.001
                                                 outfmt=6,
                                                 out=os.path.join(out_folder, "blast.tmp"),
                                                 max_target_seqs=1)
            
        else:
            blast_cline = NcbiblastnCommandline(query=blast,
                                                db=reference,
                                                evalue=0.001,
                                                outfmt=6,
                                                out=os.path.join(out_folder, "blast.tmp"))
        logging.info(f"Blast:\t{blast_cline}")
        stdout, stderr = blast_cline()

        # a,b,c = shell_command.shell_command('tblastn -query %s -db %s -evalue 1e-5 -max_target_seqs 1 -outfmt 6 > blast.tmp' % (blast, reference))
        # a,b,c = shell_command.shell_command('tblastn -query %s -db %s -evalue 1e-5 -max_target_seqs 1 -outfmt 6' % (blast, reference))

        blast2data, queries = blast_utils.remove_blast_redundancy([os.path.join(out_folder, "blast.tmp")], check_overlap=False, out_folder=self.circos_data_dir)

        o = open(os.path.join(out_folder, 'circos_blast.txt'), "w")
        l = open(os.path.join(out_folder, 'circos_blast_labels.txt'), "w")

        for contig in blast2data:
            cname = re.sub("\|", "", contig)
            for gene in blast2data[contig]:
                if float(blast2data[contig][gene][0]) >= identity_cutoff:  # 80,20
                    location = sorted(blast2data[contig][gene][1:3])
                    o.write("%s\t%s\t%s\n" % (
                    contig, location[0] + self.contigs_add[cname][0], location[1] + self.contigs_add[cname][0]))
                    l.write("%s\t%s\t%s\t%s\n" % (
                    contig, location[0] + self.contigs_add[cname][0], location[1] + self.contigs_add[cname][0], gene))

        o.close()

    def add_multiple_genome_tracks(self, track_file_list, highlight_list=[], condensed=False):
        logging.info(f'Track file list {track_file_list}')
        import os

        # r1 doit etre plus grand que r1
        r1 = self.last_track
        # r0 = self.last_track-0.015
        n = 0
        hc = 0
        for i, orthofile in enumerate(track_file_list):
            n += 1
            if orthofile not in highlight_list:
                if not condensed:
                    r0 = r1 - 0.015
                else:
                    r0 = r1 - 0.008
                if n % 2 == 0:  # orrd-9-seq # blues
                    self.circos_reference.add_plot(orthofile, type="heatmap", r1="%sr" % r1, r0="%sr" % r0,
                                                   color="ylorrd-9-seq", fill_color="", thickness="2p", z=1, rules="",
                                                   backgrounds="", url="")
                else:
                    self.circos_reference.add_plot(orthofile, type="heatmap", r1="%sr" % r1, r0="%sr" % r0,
                                                   color="gnbu-9-seq", fill_color="", thickness="2p", z=1, rules="",
                                                   backgrounds="", url="")

                # circos.add_highlight(orthofile, fill_color="ortho3", r1="%sr" % r1, r0= "%sr" % r0, href=href)
                if not condensed:
                    r1 = r1 - 0.023  # 046
                    r0 = r0 - 0.023  # 046
                else:
                    r1 = r1 - 0.011  # 046
                    r0 = r0 - 0.011  # 046
            else:
                r0 = r1 - 0.013
                n -= 1
                if hc > 9:
                    hc = 0
                color = 'set1-9-qual-%s' % hc
                hc += 1
                self.circos_reference.add_highlight(orthofile, r1="%sr" % r1, r0="%sr" % r0, fill_color=color)
                r1 = r1 - 0.020  # 046
                # r0 = r0-0.011

        try:
            self.last_track = r0
        except:
            pass
        self.config = self.circos_reference.get_file()

        # t = open('circos.config', "w")

    def add_samtools_depth_track(self, samtools_file, lower_cutoff=50, top_cutoff=5000, r0=0.8, r1=0.7):

        rules = """
        <rules>
        <rule>
        condition          = var(value) > %s
        color              = green
        fill_color         = lgreen
        </rule>

        <rule>
        condition          = var(value) < %s
        color              = red
        fill_color         = lred
        </rule>
        </rules>

        """ % (top_cutoff,
               lower_cutoff)

        self.circos_reference.add_plot(samtools_file,
                                       type="histogram",
                                       r1="%sr" % (r1),
                                       r0="%sr" % (r0),
                                       color="black",
                                       fill_color="grey_a5",
                                       thickness="1p",
                                       z=1,
                                       rules=rules,
                                       backgrounds="",
                                       url="")

        self.config = self.circos_reference.get_file()

    def run_circos(self, config_file="circos.config", out_prefix="circos"):
        from mummer2circos import shell_command
        cmd = 'circos -outputfile %s.svg -conf %s' % (out_prefix, config_file)
        stdout_str, stderr_str, returncode = shell_command.shell_command(cmd)
        with open("circos.log", "w") as f:
            f.write(stdout_str.decode("utf-8") )
        with open("circos.err", "w") as f:
            f.write(stderr_str.decode("utf-8") )
        if returncode == 0:
            logging.info(f"Circos plot generated sucessfully, see {out_prefix}.svg & {out_prefix}.png")
            if self.heatmap:
                main_conf = ["- main configuration file:".ljust(35), "circos.config"]
                gc_skew = ["- GC skew:".ljust(35), "circos_GC_skew.txt"]
                gc_content = ["- GC content:".ljust(35), "circos_GC_content.txt"]
                contigs_def = ["- Contigs definition:".ljust(35), "circos_contigs.txt"]
                circos_files = [main_conf,contigs_def, gc_skew, gc_content]
                if self.gbk2orf:
                    circos_files.append(["- Reference genome ORFs strand -:".ljust(35), "circos_orf_minus.txt"])
                    circos_files.append(["- Reference genome ORFs strand +:".ljust(35), "circos_orf_plus.txt"])
                for n, heatmap_track in enumerate(self.genome_list):
                    circos_files.append([f"- Heatmap track {n+1}:".ljust(35), heatmap_track])
                if self.blast:
                    circos_files.append([f"- Blast hits labels:".ljust(35), "circos_blast_labels.txt"])
                if self.samtools_depth:
                    for n, depth_track in enumerate(self.samtools_depth_track_list):
                        circos_files.append([f"- Depth track {n+1}:".ljust(35), depth_track])
                log_str = ['\t'.join(i) for i in circos_files]
                logging.info(f"Circos files saved in {self.circos_data_dir}: \n\t" + '\n\t'.join(log_str))
            logging.info(f'The plot can be modified (to change color scales, remove or add tracks...) by modifining the file "circos.config" and excuting the command "circos -conf circos.config"')
        else:
            raise CircosException("Circos problem, check circos log file (circos.log and circos.err)... quitting")

    def clean_tmp_files(self):
        import os
        d = os.getcwd()
        purge(d, ".*.ndb")
        purge(d, ".*.nhr")
        purge(d, ".*.nin")
        purge(d, ".*.njs")
        purge(d, ".*.not")
        purge(d, ".*.nsq")
        purge(d, ".*.ntf")
        purge(d, ".*.nto")

    def write_circos_files(self, config, color_config):
        import os
        with open(f"{self.circos_data_dir}/circos.config", 'w') as f:
            f.write(config)
        ''''
        try:
            os.mkdir(os.path.join(self.working_dir, 'etc/'))
        except:
            pass
        with open("etc/brewer.all.conf", "w") as f:
            f.write(color_config)
        '''

    def id_generator(self, size=6, chars=string.ascii_lowercase):  # + string.digits
        return ''.join(random.choice(chars) for _ in range(size))

    def get_circos_config(self, last_seq_id, first_seq_id, mid1, mid2, link=True, heat=False):

        link_code = '''

        <links>

        <link>
        ribbon = yes
        file          = circos.link
        color         = orrd-9-seq
        radius        = 0.95r
        bezier_radius = 0.1r
        thickness     = 3

        #<rules>
        #<rule>
        # you can also test whether only one end is
        # reversed using var(inv)
        #condition  = var(rev1) && ! var(rev2)
        #color      = blue
        #</rule>
        #<rule>
        #condition  = var(rev2) && ! var(rev1)
        #color      = blue
        #</rule>
        #<rule>
        #condition  = var(rev1) && var(rev2)
        #color      = orange
        #</rule>

        #</rules>

        </link>

        </links>


        '''

        plot_templat = '''

        <plots>

        %s

        </plots>



        '''

        chr_spacing = '''

           <pairwise %s %s>
         spacing = 14u
        </pairwise>

        '''

        circos_config = f'''

         karyotype = {self.circos_data_dir}/circos_contigs.txt
         chromosomes_units           = 10000
         chromosomes_display_default = yes
        <ideogram>

         <spacing>
         default            = %su

         %s

         </spacing>



         # thickness and color of ideograms
         thickness          = 20p
         stroke_thickness   = 1
         stroke_color       = black

         # the default chromosome color is set here and any value
         # defined in the karyotype file overrides it
         fill               = yes
         fill_color         = black

         # fractional radius position of chromosome ideogram within image
         radius             = 0.85r
         show_label         = no
         label_font         = default
         label_radius       = dims(ideogram,radius) + 0.175r
         label_size         = 30
         label_parallel     = no

         # show_bands determines whether the outline of cytogenetic bands
         # will be seen
         show_bands         = yes
         band_stroke_thickness = 1
          # in order to fill the bands with the color defined in the karyotype
         # file you must set fill_bands
         fill_bands         = yes
         band_transparency  = 1

         </ideogram>
        show_ticks         = yes
         show_tick_labels   = yes

         <ticks>

         tick_label_font    = condensed
         radius             = dims(ideogram,radius_outer)
         label_offset       = 8p
         label_size         = 8p
         color              = black
         thickness          = 4p

         #<tick>
         #spacing           = 100u
         #size              = 16p
         #label_multiplier  = 1e-3
         #show_label        = yes
         #label_size        = 35p
         #format            = %%d kb
         #thickness         = 5p
         #</tick>


         <tick>
         skip_first_label = yes
         multiplier   = 1/1u
         spacing           = 1u
         size              = 15p
         show_label        = no
         label_size        = 25p
         </tick>

         <tick>
         skip_first_label = yes
         multiplier   = 10/1u
         spacing           = 10u
         size              = 10p
         show_label        = yes
         label_size        = 25p
         format            = %%.1d kb
         </tick>

        <tick>
        multiplier   = 1
        position = end
        size              = 16p
        show_label        = yes
        label_size        = 25p
        format    = %%d bp
        </tick>

        <tick>
        position = 0u
        size              = 16p
        label_size        = 25p
        label = a
        </tick>

         </ticks>


        %s



        <colors>
         <<include colors.my>>
         #<<include brewer.all.conf>>
         </colors>
         <image>
         image_map_use      = yes
         image_map_overlay  = yes
         image_map_overlay_stroke_color     = red
         <<include etc/image.conf>>
         </image>
         <<include etc/colors_fonts_patterns.conf>>
         # system and debug settings
         <<include etc/housekeeping.conf>>
         anti_aliasing*     = no

         '''

        if link:

            if last_seq_id == mid1 and mid2 == first_seq_id:
                ch = chr_spacing % (last_seq_id, first_seq_id)
                return circos_config % (0.5, ch, link_code)
            else:
                ch1 = chr_spacing % (last_seq_id, first_seq_id)
                ch2 = chr_spacing % (mid1, mid2)
                return circos_config % (0.5, ch1 + ch2, link_code)

    def get_karyotype_from_fasta(self,
                                 fasta1,
                                 fasta2,
                                 hit_list,
                                 query_list,
                                 filter_ref=True,
                                 filter_query=True,
                                 out="circos_contigs.txt",
                                 both_fasta=True,
                                 cumul=True):
        from Bio import SeqIO
        import re

        fasta_data1 = [i for i in SeqIO.parse(open(fasta1), "fasta")]
        if fasta2:
            fasta_data2 = [i for i in SeqIO.parse(open(fasta2[0]), "fasta")]

        self.contig2start_psoition = {}

        with open(out, 'w') as f:
            # chr - Rhab Rhab 0 1879212 spectral-5-div-4
            i = 0
            contig_start = 0
            contig_end = 0
            for record in fasta_data1:
                name = re.sub("\|", "", record.id)
                # cumulated length if not link plot and not filter_ref (TODO, put it as an argument?)
                if cumul:
                    contig_start = contig_end + 1
                contig_end = contig_start + len(record)

                # keep in memory
                self.contig2start_psoition[name] = contig_end + 1

                if i == 4:
                    i = 0
                i += 1

                if filter_ref:
                    if name in hit_list:
                        n4 = name
                        if not 'n2' in locals():
                            n2 = name
                        # spectral-5-div-%s

                        f.write("chr - %s %s %s %s greys-3-seq-%s\n" % (name, name, contig_start, contig_end, i))  # i,
                else:
                    n4 = name
                    if not 'n2' in locals():
                        n2 = name
                    # spectral-5-div-%s

                    f.write("chr - %s %s %s %s greys-3-seq-%s\n" % (name, name, contig_start, contig_end, i))  # , i

            if both_fasta:
                for record in fasta_data2:
                    name = re.sub("\|", "", record.id)
                    if filter_query:
                        if name in query_list:
                            f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (name, name, 0, len(record), i))
                    else:
                        f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (name, name, 0, len(record), i))
        if fasta2:
            n1 = re.sub("\|", "", fasta_data2[-1].name)
            n3 = re.sub("\|", "", fasta_data2[0].name)
            if not 'n2' in locals():
                n2 = n1
            return (n1, n2, n3, n4)
        else:
            return None

    def execute_megablast(self, fasta1, fasta2):
        import os
        from mummer2circos import shell_command
        for one_fasta in fasta2:
            out_prefix = os.path.basename(one_fasta).split('.')[0]

            cmd1 = "makeblastdb -in %s -dbtype nucl" % one_fasta
            cmd2 = 'blastn -task megablast -query %s -db %s -evalue 1e-5 -outfmt 6 -out circos_tmp/blast_result_%s.tab' % (fasta1,
                                                                                                                one_fasta,
                                                                                                                out_prefix)
            a, b, c = shell_command.shell_command(cmd1)
            a, b, c = shell_command.shell_command(cmd2)
            if c != 0:
                raise Exception(f"Problem with blast:\n{a}\n{b}")

    def justLinks(self, coords_input):
        ##Find first row after header
        i = 0
        for row in coords_input:
            if len(re.split(r'\t+', row)) > 3:
                headerRow = i + 1
                break
            else:
                i += 1

        return coords_input[headerRow:None]

    def get_link(self, coords_input, link_file="circos.link", algo="nucmer"):
        import re
        import matplotlib.cm as cm
        import matplotlib as mpl

        if algo == 'promer':
            shift = 4
        elif algo == 'nucmer':
            shift = 0

        with open(coords_input, 'rU') as infile:
            rawLinks = self.justLinks(infile.readlines());

        all_id = [float(re.split(r'\t+', i.rstrip('\n'))[6]) for i in rawLinks]
        all_start = []
        all_ends = []

        id_min = min(all_id)
        id_max = max(all_id)

        norm = mpl.colors.Normalize(vmin=id_min, vmax=id_max)
        cmap = cm.Blues
        cmap_blue = cm.Reds  # OrRd

        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
        c = open('colors.my', 'w')

        with open(link_file, 'w') as f:
            # c.write()
            i = 1
            hit_list = []
            query_list = []
            for row in rawLinks:
                color = ''
                l = re.split(r'\t+', row.rstrip('\n'))

                if int(l[0]) < int(l[1]) and int(l[2]) > int(l[3]):
                    color = m.to_rgba(float(l[6]))

                elif int(l[0]) > int(l[1]) and int(l[2]) < int(l[3]):
                    color = m.to_rgba(float(l[6]))

                else:
                    color = m2.to_rgba(float(l[6]))

                color_id = self.id_generator()

                c.write('%s = %s,%s,%s,%s\n' % (color_id,
                                                int(round(color[0] * 250, 0)),
                                                int(round(color[1] * 250, 0)),
                                                int(round(color[2] * 250, 0)),
                                                0.5))

                f.write(re.sub("\|", "", l[9 + shift]) + '\t' + l[0] + '\t' + l[1] + '\t' + re.sub("\|", "", l[
                    10 + shift]) + '\t' + l[2] + '\t' + l[3] + '\tcolor=%s' % color_id + '\n')
                i += 1
                hit_list.append(re.sub("\|", "", l[9 + shift]))
                query_list.append(re.sub("\|", "", l[10 + shift]))
        return (hit_list, query_list)

    def megablast2heatmap(self, megablast_input, link_file="circos.heat", col=250):
        import re
        import os
        import matplotlib.cm as cm
        import matplotlib as mpl

        with open(megablast_input, 'rU') as infile:
            rawLinks = [i.rstrip().split('\t') for i in infile]

        all_id = [float(i[2]) for i in rawLinks]

        try:
            id_min = min(all_id)
        except:
            return None
        id_max = max(all_id)

        norm = mpl.colors.Normalize(vmin=id_min, vmax=id_max)
        cmap = cm.Blues
        cmap_blue = cm.Reds  # OrRd
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
        c = open('colors.my', 'w')

        contig2start_stop_list = {}

        heatmap_file = os.path.basename(megablast_input).split('.')[0] + '.heat' 
        heatmap_file = re.sub("blast_result_", "", heatmap_file)
        with open(heatmap_file, 'w') as f:
            i = 1
            hit_list = []
            query_list = []
            for row in rawLinks:
                color = ''
                if row[0] not in contig2start_stop_list:
                    contig2start_stop_list[row[0]] = {}
                    contig2start_stop_list[row[0]]["start"] = [row[6]]
                    contig2start_stop_list[row[0]]["stop"] = [row[7]]
                else:
                    contig2start_stop_list[row[0]]["start"].append(row[6])
                    contig2start_stop_list[row[0]]["stop"].append(row[7])

                color = m2.to_rgba(float(row[2]))

                color_id = self.id_generator()


                c.write('%s = %s,%s,%s,%s\n' % (color_id,
                                                int(round(color[0] * col, 0)),
                                                int(round(color[1] * col, 0)),
                                                int(round(color[2] * col, 0)),
                                                0.5))


                f.write(re.sub("\|", "", row[0]) + '\t' + row[6] + '\t' + row[7] + '\t' + row[2] + "\tz=%s\t" % row[
                    2] + '\n')

                i += 1
                hit_list.append(re.sub("\|", "", row[0]))
                query_list.append(re.sub("\|", "", row[1]))
        return (hit_list, query_list, contig2start_stop_list)

    def gbk2circos_data(self, gbk_file,
                        minus_file="circos_orf_minus.txt",
                        plus_file="circos_orf_plus.txt"):

        m = open(minus_file, 'w')
        p = open(plus_file, 'w')
        from Bio import SeqIO

        with open(gbk_file, 'r') as f:
            for record in SeqIO.parse(f, 'genbank'):
                for feature in record.features:
                    start = str(feature.location.start + self.contigs_add[record.id][0])
                    end = str(feature.location.end + self.contigs_add[record.id][0])
                    if abs(int(end) - int(start)) > 100000:
                        continue
                    if feature.type == 'CDS':
                        # deal with impossible size ORF
                        # if 'hyp' in feature.qualifiers['product'][0]:
                        #    col='blue'
                        # else:
                        col = 'grey_a1'

                        if self.locus2taxo_color:
                            try:
                                col = self.locus2taxo_color[feature.qualifiers['locus_tag'][0]]
                            except KeyError:
                                raise('Locus not found. Taxonomy file shoud report locus_tag2taxonomy correspondance.')

                        if int(feature.location.strand) == 1:

                            p.write('%s\t%s\t%s\tfill_color=%s\n' % (record.id,
                                                                     start,
                                                                     end,
                                                                     col))
                        else:
                            m.write('%s\t%s\t%s\tfill_color=%s\n' % (record.id,
                                                                     start,
                                                                     end,
                                                                     col))

                    elif feature.type == 'rRNA' or feature.type == 'tRNA' or feature.type == 'tmRNA' or feature.type == 'ncRNA':
                        if feature.location.strand == 1:
                            p.write('%s\t%s\t%s\tfill_color=red\n' % (record.id,
                                                                      start,
                                                                      end))
                        else:
                            m.write('%s\t%s\t%s\tfill_color=red\n' % (record.id,
                                                                      start,
                                                                      end))
                    else:
                        pass
        m.close()
        p.close()

        return minus_file, plus_file

    def samtools_depth2circos_data(self, samtools_depth_file, i, window=1000):
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
        with open(f'{self.circos_data_dir}/circos_samtools_depth_%s.txt' % i, 'w') as g:
            for contig in contig2coverage:
                # split list by chunks of 1000
                mychunks = [i for i in chunks(contig2coverage[contig], window)]
                for i, cov_list in enumerate(mychunks):
                    median_depth = numpy.median(cov_list)
                    if median_depth > (2 * all_contigs_median):
                        median_depth = (2 * all_contigs_median) + 1
                    try:
                        g.write("%s\t%s\t%s\t%s\n" % (contig, (i * window) + self.contigs_add[contig][0],
                                                  ((i * window) + window - 1) + self.contigs_add[contig][0],
                                                  median_depth))
                    except:
                        # contig present in depth file missing from fasta file
                        # skip as the case might happen (eg mapping before contig filtering based on depth)
                        continue
        return all_contigs_median

