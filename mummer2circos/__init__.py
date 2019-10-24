
def main():
    import mummer2circos
    import argparse
    ###Argument handling.
    arg_parser = argparse.ArgumentParser(description='');
    # arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-r", "--fasta1", help="reference fasta")
    arg_parser.add_argument("-q", "--fasta2", help="query fasta", nargs='+')
    arg_parser.add_argument("-fr", "--filterr", action="store_true",
                            help="do not remove reference sequences without any similarity from the plot (default False)")
    arg_parser.add_argument("-fq", "--filterq", action="store_false",
                            help="do not remove query sequences without any similarity from the plot (default False)")
    arg_parser.add_argument("-l", "--link", help="link circos and not heatmap circos", action="store_true")
    arg_parser.add_argument("-s", "--samtools_depth", help="samtools depth file", nargs="+")
    arg_parser.add_argument("-o", "--output_name", help="output circos pefix", default="nucmer2circos")
    arg_parser.add_argument("-g", "--gaps", help="highlight gaps", action="store_true")
    arg_parser.add_argument("-gb", "--genbank", help="add ORF based on GBK file", default=False)
    arg_parser.add_argument("-b", "--blast", help="highlight blast hits (-outfmt 6)")
    arg_parser.add_argument("-bc", '--blast_identity_cutoff', type=int, help="Blast identity cutoff", default=80)
    arg_parser.add_argument("-n", "--highlight", help="highlight instead of heatmap corresponding list of records",
                            nargs="+")
    arg_parser.add_argument("-a", "--algo", help="algorythm to use to compare the genome (megablast, nucmer or promer)",
                            default="nucmer")
    arg_parser.add_argument("-m", "--min_gap_size", help="minimum gap size to consider", default=1000)
    arg_parser.add_argument("-bn", '--blastn', action="store_true", help="excute blastn and not blastp")
    arg_parser.add_argument("-w", '--window', type=int, help="window size (default=1000)", default=1000)
    arg_parser.add_argument("-ss", '--secretion_systems', type=str, help="macsyfinder table", default=False)
    arg_parser.add_argument("-c", '--condensed', action="store_true", help="condensed display (for mor tracks)",
                            default=False)
    arg_parser.add_argument("-lf", '--label_file', type=str,
                            help="label file ==> tab file with: contig, start, end label (and color)", default=False)
    arg_parser.add_argument("-lt", '--locus_taxonomy', type=str,
                            help="Color locus based on taxonomy: tab delimited file with: locus\tphylum. " \
                                 " Color locus matching the Taxon set in comment as the first row (#Chlamydiae)", default=False)


    args = arg_parser.parse_args()

    if args.highlight is None:
        args.highlight = []
    ###Variable Definitions

    ##Run main
    if args.fasta2 is None:
        fasta2=[]
    else:
        fasta2=args.fasta2
    circosf = mummer2circos.Fasta2circos(args.fasta1,
                           fasta2,
                           args.filterr,
                           args.filterq,
                           heatmap=args.link,
                           samtools_depth=args.samtools_depth,
                           gaps=args.gaps,
                           blast=args.blast,
                           highlight_list=args.highlight,
                           algo=args.algo,
                           min_gap_size=int(args.min_gap_size),
                           blastn=args.blastn,
                           gbk2orf=args.genbank,
                           window_size=args.window,
                           secretion_systems=args.secretion_systems,
                           condensed_tracks=args.condensed,
                           label_file=args.label_file,
                           locus_taxonomy=args.locus_taxonomy,
                           blast_identity_cutoff=args.blast_identity_cutoff)

    circosf.write_circos_files(circosf.config, circosf.brewer_conf)
    circosf.run_circos(out_prefix=args.output_name)


