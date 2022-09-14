
import logging


def run_prodigal(fasta_seq, output_name='temp.faa'):
    from Bio import SeqIO
    from mummer2circos import shell_command
    from io import StringIO
    from tempfile import NamedTemporaryFile
    # -q quiet
    # -a Write protein translations to the selected file
    # -i Specify input file
    # -c:  Closed ends.  Do not allow genes to run off edges. # not activated
    cmd = "prodigal -q -a %s -i %s" % (output_name, fasta_seq)

    sdt_out, sdt_err, err = shell_command.shell_command(cmd)
    logging.info (sdt_out)
    logging.info  (sdt_err)
    stdout_str, stderr_str, err_code = shell_command.shell_command('sed -i "s/*//g" %s' % output_name)
    
    if err_code != 0:
        raise IOError("Problem formating BLAST database:", stderr_str)
    
    return output_name


class Hmm():

    def __init__(self, hmm_profiles, database, output_dir=False, call_genes=False, score_cutoff=100, cut_tc=False, filter_bitscore='best'):
        from tempfile import NamedTemporaryFile

        assert isinstance(hmm_profiles, list)

        if isinstance(database, list) and len(database)>1:
            # if multiple databases, compare the bitscore obtained for each profile against each protein
            # use either the best biscrore or the median bitscrore to filter results
            self.bitscore_litering = filter_bitscore
            self.multiple_databases = True
            self.profile2scores = {}

        self.hmm_profiles = hmm_profiles
        self.database = database
        # -T <x>     : report sequences >= this score threshold in output
        # out, query and db
        #self.hmmer_cmd = 'hmmsearch -T %s -E 1e-10 -o %s %s %s'
        if not cut_tc:
            self.hmmer_cmd = 'hmmsearch -T %s -o' % score_cutoff + ' %s %s %s'
        else:
            self.hmmer_cmd = 'hmmsearch --cut_tc -o %s %s %s'
        self.hmmer_score_cutoff = score_cutoff
        self.hmmer_output_list = []

        if call_genes:
            #if not output_dir:

            temp_prodigal_file = NamedTemporaryFile()
            self.database = temp_prodigal_file.name
            # add content to temporary file
            #temp_file.write(str(self))

            run_prodigal(database, output_name=self.database)

    def run_hmmer(self, profiles=False):
        from tempfile import NamedTemporaryFile
        from mummer2circos import shell_command

        if not profiles:
            profiles = self.hmm_profiles

        header = ["profile_id",
                "profile_length",
                "best_hit_id",
                "bias",
                "bitscore",
                "evalue",
                "query_start",
                "query_end",
                "query_coverage",
                "hit_start",
                "hit_end"]
        results = []#[header]
        for profile in profiles:
            temp_file = NamedTemporaryFile()
            self.hmmer_output_list.append(temp_file.name)
            if not isinstance(self.database, list):
                cmd = self.hmmer_cmd % (temp_file.name, profile, self.database)

                stout, sterr, code = shell_command.shell_command(cmd) # self.hmmer_score_cutoff,
                if code != 0:
                    raise IOError("Error", sterr)

                parsed_data = self._parse_hmmsearch(temp_file.name)

                if len(parsed_data) == 0:
                    logging.info ('No domains respecting score threshold for %s, continue...' % profile)
                    continue

                if not isinstance(parsed_data[0], dict):
                    results.append(['%s' % parsed_data[0], '-', '-', '-', '-', '-', '-', '-', '-', '-'])
                else:
                    hsp_list = parsed_data
                    for x in range(0,len(hsp_list)):
                        #results += '\t'.join([str(hsp_list[x][i]) for i in header])
                        #results += '\n'
                        results.append([str(hsp_list[x][i]) for i in header])
            else:
                # multiple databases: performing bitscore filtering
                self.biodb2best_hits = {}
                for database in self.database:
                    stout, sterr, code = shell_command.shell_command(self.hmmer_cmd % (self.hmmer_score_cutoff, temp_file.name, profile, self.database))
                    if code != 0:
                        raise IOError("Error", sterr)
                        #import sys
                        #sys.stdout.write("\n%s\n%s\n" % (stout, sterr))
                        #sys.exit()

                    parsed_data = self._parse_hmmsearch(temp_file.name)


                    '''
                    if not isinstance(parsed_data[0], dict):
                        pass
                    else:
                        # all hsp have the same bitscore, only use the first hsp
                        if parsed_data[0]['profile_id'] not in self.profile2scores:
                            self.profile2scores[parsed_data[0]['profile_id']] = [parsed_data[0]['bitscore']]
                        else:
                            self.profile2scores[parsed_data[0]['profile_id']].append(parsed_data[0]['bitscore'])
                        hsp_list = parsed_data
                        for x in range(0,len(hsp_list)):
                            results += '\t'.join([str(hsp_list[x][i]) for i in header])
                            results += '\n'
                    '''

        return results

    def _parse_hmmsearch(self, hmmsearch_result):
        from Bio import SearchIO

        result_handle = open(hmmsearch_result, 'r')
        hmmer_records = [i for i in SearchIO.parse(result_handle, 'hmmer3-text')]

        try:
            best_hit_id = hmmer_records[0].hits[0].id
        except IndexError:
            return [hmmer_records[0].id]
        else:

            hsp_list = []
            profile_id =  hmmer_records[0].id
            '''
            'append', 'bias', 'bitscore', 'description', 'description_all', 'domain_exp_num', 'domain_obs_num',
            'evalue', 'filter', 'fragments', 'hsps', 'id',
            'id_all', 'index', 'is_included', 'map', 'pop', 'query_description', 'query_id', 'sort'

            '''
            for hsp in hmmer_records[0].hits[0].hsps:
                result = {
                    "profile_id" : hmmer_records[0].id,
                    "profile_length" : hmmer_records[0].seq_len,
                    "best_hit_id" : hmmer_records[0].hits[0].id,
                    "bias" : hmmer_records[0].hits[0].bias,
                    "bitscore" : hmmer_records[0].hits[0].bitscore,
                    "hit_start" : hsp.hit_start,
                    "hit_end" : hsp.hit_end,
                    "query_start" : hsp.query_start,
                    "query_end" : hsp.query_end,
                    "evalue" : hmmer_records[0].hits[0].evalue,
                    "query_coverage" : round((hsp.query_end-hsp.query_start)/float(hmmer_records[0].seq_len),2)
                    }
                hsp_list.append(result)
            return hsp_list

    def filter_hmmer_results(self):
        pass


def remove_blast_redundancy(blast_file_list, 
                            check_overlap=True, 
                            out_folder='.'):

    import os
    '''

    0   Query
    1   Subject
    2   % id
    3   alignment length
    4   mistmatches
    5   gap openings
    6   q.start
    7   q.end
    8   s.start
    9  s.end
    10  e-value
    11  bit score

    :param blast_file_list:
    :param check_overlap:
    :return:
    '''

    blast2data = {}
    queries = []

    # keep only best hit
    # todo make it proper
    for one_blast_file in blast_file_list:
        result_handle = open(one_blast_file, 'r')
        best_hit_handle = open(os.path.join(out_folder, "niq_best.tmp"), 'w')
        hit_list = []
        for line in result_handle:
            if line.split('\t')[0] in hit_list:
                continue
            else:
                hit_list.append(line.split('\t')[0])
                best_hit_handle.write(line)
        best_hit_handle.close()

        with open(os.path.join(out_folder, "niq_best.tmp"), 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[1] not in blast2data:
                    blast2data[line[1]] = {}
                    blast2data[line[1]][line[0]] = [float(line[2]), int(line[8]), int(line[9])]
                else:
                     blast2data[line[1]][line[0]] = [float(line[2]), int(line[8]), int(line[9])]
                if line[0] not in queries:
                    queries.append(line[0])

    if check_overlap:
        for one_blast in blast2data.keys():
            for ref_gene in blast2data[one_blast].keys():

                for query_gene in blast2data[one_blast].keys():
                    overlap = False
                    if ref_gene == query_gene:
                        continue
                    # check if position is overlapping
                    try:
                        sorted_coordinates = sorted(blast2data[one_blast][ref_gene][1:3])
                        if blast2data[one_blast][query_gene][1] <= sorted_coordinates[1] and blast2data[one_blast][query_gene][1]>= sorted_coordinates[0]:
                            overlap =True
                        sorted_coordinates = sorted(blast2data[one_blast][query_gene][1:3])
                        if blast2data[one_blast][ref_gene][1] <= sorted_coordinates[1] and blast2data[one_blast][ref_gene][1]>= sorted_coordinates[0]:
                            overlap =True
                        if overlap:
                            if float(blast2data[one_blast][ref_gene][0]) > float(blast2data[one_blast][query_gene][0]):
                                del blast2data[one_blast][query_gene]
                            else:
                                del blast2data[one_blast][ref_gene]
                                break
                    except KeyError:
                        # colocation already resolved
                        pass
                         
    return blast2data, queries



class Blast():

    def __init__(self, query, database, protein=False, formatdb=False, best_hit_only=True):
        import os
        from Bio import SeqRecord, SeqIO
        
        if type(query) == list or isinstance(query, SeqRecord.SeqRecord):
            from io import StringIO
            from tempfile import NamedTemporaryFile
            temp_query = NamedTemporaryFile(delete=False, mode='w')
            fastastr = StringIO()
            SeqIO.write(query, fastastr, 'fasta')
            temp_query.write(fastastr.getvalue())
            temp_query.flush()
            self.query = temp_query.name
            # add content to temporary file

        elif type(query) == str:
            self.query = query
            self.query = query

        else:
            raise TypeError('wrong inut format: either SeqRecord or string')


        if type(database) == list or isinstance(database, SeqRecord.SeqRecord):
            from io import StringIO
            from tempfile import NamedTemporaryFile
            temp_db = NamedTemporaryFile(delete=False, mode='w')
            fastastrdb = StringIO()
            SeqIO.write(database, fastastrdb, 'fasta')
            temp_db.write(fastastrdb.getvalue())
            temp_db.flush()
            self.database = temp_db.name
            # add content to temporary file

        elif type(database) == str or type(database) == unicode:
            self.database = database

        else:
            raise TypeError('wrong inut format: either SeqRecord or string')

        self.protein = protein
        self.best_hit_only = best_hit_only
        self.formatdb = formatdb
        self.working_dir = os.getcwd()
        self.blast_path_var= 'export BLASTDB=/tmp/'
        import shell_command
        stdout_str, stderr_str, err_code = shell_command.shell_command(self.blast_path_var)
        if err_code != 0:
            raise IOError("Error", stderr_str)


    def id_generator(self, size=6, chars=False): # + string.digits
        import random
        import string

        if not chars:
            chars = string.ascii_lowercase
        return ''.join(random.choice(chars) for _ in range(size))

    def format_database(self):
        from mummer2circos import shell_command

        new_database = self.id_generator(8)

        if self.protein:
            cmd = 'makeblastdb -in %s -dbtype prot' % (self.database)
            stdout_str, stderr_str, err_code = shell_command.shell_command(cmd )
        else:
            cmd = 'makeblastdb -in %s -dbtype nucl' % (self.database)
            stdout_str, stderr_str, err_code = shell_command.shell_command(cmd)

        if err_code != 0:
            raise IOError("Problem formating BLAST database:", stderr_str)

        self.database_path = '/tmp/%s.temp' % new_database


    def run_blastp(self):
        from Bio.Blast.Applications import NcbiblastpCommandline
        import os

        blast_id = self.id_generator(8)

        outpath = os.path.join('/tmp/%s.tab' % blast_id)
        
        blastp_cline = NcbiblastpCommandline(query= self.query,
                                            db=self.database,
                                            evalue=0.005,
                                            outfmt=6,
                                            out=outpath,
                                            num_threads=8)
        logging.info (blastp_cline)
        stdout, stderr = blastp_cline()
        logging.info (stderr)

        with open(outpath, 'r') as result_handle:

            self.best_hit_list = []
            self.complete_hit_list = []
            for line in result_handle:
                self.complete_hit_list.append(line.rstrip().split('\t'))
                if line.split('\t')[0] in self.best_hit_list:
                    continue
                else:
                    self.best_hit_list.append(line.rstrip().split('\t'))

        return outpath


    def run_blastn(self, min_identity=50):
        from Bio.Blast.Applications import NcbiblastnCommandline
        import os

        blast_id = self.id_generator(8)

        outpath = os.path.join('/tmp/%s.tab' % blast_id)
        
        blastn_cline = NcbiblastnCommandline(query= self.query,
                                             task="dc-megablast",
                                            db=self.database,
                                            evalue=0.1,
                                            outfmt=6,
                                            out=outpath,
                                            num_threads=8,
                                            max_hsps=1000,
                                             perc_identity=min_identity)

        logging.info(blastn_cline)
        stdout, stderr = blastn_cline()
        logging.info(stderr)

        with open(outpath, 'r') as result_handle:

            self.best_hit_list = []
            self.complete_hit_list = []
            for line in result_handle:
                self.complete_hit_list.append(line.rstrip().split('\t'))
                if line.split('\t')[0] in self.best_hit_list:
                    continue
                else:
                    self.best_hit_list.append(line.rstrip().split('\t'))

        return outpath

    def run_tblastx(self,evalue=0.1):
        from Bio.Blast.Applications import NcbitblastxCommandline
        import os

        blast_id = self.id_generator(8)

        outpath = os.path.join('/tmp/%s.tab' % blast_id)
        
        blastn_cline = NcbitblastxCommandline(query= self.query,
                                            db=self.database,
                                            evalue=evalue,
                                            outfmt=6,
                                            out=outpath,
                                            num_threads=8,
                                            max_hsps=1000)

        logging.info(blastn_cline)
        stdout, stderr = blastn_cline()
        logging.info(stderr)

        with open(outpath, 'r') as result_handle:

            self.best_hit_list = []
            self.complete_hit_list = []
            for line in result_handle:
                self.complete_hit_list.append(line.rstrip().split('\t'))
                if line.split('\t')[0] in self.best_hit_list:
                    continue
                else:
                    self.best_hit_list.append(line.rstrip().split('\t'))

        return outpath

    def run_tblastn(self):


        '''
        tblastn_cline = NcbitblastnCommandline(query='dnaa.temp',
                                             db=contig_file,
                                             evalue=0.001,
                                             outfmt=5,
                                             out="dnaa_blast2.xml")

        stdout, stderr = tblastn_cline()

        result_handle = open("dnaa_blast2.xml", 'r')
        blast_records = [i for i in NCBIXML.parse(result_handle)]
        '''
        pass
