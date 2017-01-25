import argparse, sys, os
from queryUtils import GEO_query
from GEOsearch import SOFTQuickParser
from Related_Sample_Search import Related_Sample_Search
from setup import get_settings, setup
from update import update


def Help():
    print "\nGEO Chip Finder"
    print "A list of functions for GEO Chip Seq Sample Finder, please try:\npython GCF.py -h"
    print "\nFuctions:"
    print "\tsearch:\n\tsearch chip-seq samples and corresponding input from GEO based on the key words."
    print "\tmatch:\n\tmatch different types of samples from GEO based, for example, looking for the corresponding H3K27me3 samples for each H3K4me3 samples."
    print "\tquery:\n\tget SRR sequencing information by several common identifiers from GEO such as GSE, GSM, SRR, SRX, SRP, etc"


def GCF_search():
    '''
    this function provie an entrance to search function

    '''
    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\npython GCF.py search [optional arguments] <features>\n\nfor more help, please try: python GCF.py search -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython GCF.py search [optional arguments] <features>\n\n",
                                     description='',epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'search' to looking for samples with key words")

    parser.add_argument('feature_key_words', default=None,
                        help="list of feature key words need to used to looking for the NGS sequencing samples, "
                             "different key words need to be separated by ','")
    parser.add_argument('output_prefix', default=None, help="specify the output file prefix.")
    parser.add_argument('output_path', default=None, help="specify the output file location.")

    ## optional parameters
    parser.add_argument('-b', dest='keywords_begin', metavar='', default='',
                        help="list of key words in features need to be used to occur in the beginning of a word, "
                             "different key words need to be separated by ','")
    parser.add_argument('-t', '--type', dest='type_seq', metavar='', default='chip-seq', help="type of sequencing specified in the search. Default is 'chip-seq")
    parser.add_argument('-c', '--ignorecase', dest='ignorecase', metavar='', default=1, type=int,
                        help="specify whether case (A vs a) need to be ignored in the search. Default is 1 which means case will be ignored in the search. Set to 0 if don't want to ignore the case. ")
    parser.add_argument('--geo', dest='geo', default=0, type=int, metavar='',
                        help="specify whether search will be limited to GEO website search result. Default is 0. Set to 1 if want to perform all GEO data search.")
    parser.add_argument('--candidateslist', dest='geo_file', default=None, metavar='',
                        help="specify the file path of GSM ID list if '--hasCandidates' set to 1")
    parser.add_argument('-s','--species', dest='species', default='Homo sapiens', metavar='',
                        help="specify the samples' species. Default is Homo sapiens. Please use the species official name. For example, human is Homo sapiens."
                             "If the species name contains space, surround the name with double quotes, for example \"Homo sapiens\"")
    parser.add_argument('-e', '--encode', dest='encode_remove', default=1, type=int, metavar='',
                        help="specify whether need to remove Encode data. Default is 1. Set to 0 to keep Encode data from search.")
    parser.add_argument('-r', '--roadmap', dest='roadmap_remove', default=1, type=int, metavar='',
                        help="specify whether need to remove Roadmap data. Default is 1. Set to 0 to keep Roadmap data from search.")
    parser.add_argument('-m', '--metadata', dest='MetaData', default=None, metavar='',
                        help="specify the GSMs metadata files path")
    parser.add_argument('-p', '--process', dest='process', default=20, type=int, metavar='',
                        help="specify the number of parallel search processes want to use.")

    args = None

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print "\GEO chip seq finder\n"
        parser.print_help()
        print "\n"
        return 0
    elif len(sys.argv) >= 3:
        try:
            args = parser.parse_args()
        except:
            print "\nfor more help, please try: python GCF.py search -h\n"
            return 1

    if args is not None:
        if os.path.exists(args.output_path):
            pass
        else:
            print "Output path is not exist"
            return 1

        settings = get_settings()
        encode_pkl = settings['Encode']
        roadmap_pkl = settings['Roadmap']
        GSMGSE_pkl = settings['GSMGSE_pkl_path']

        keywords = args.feature_key_words.split(",")
        output_prefix = args.output_prefix
        output_path = args.output_path

        if args.keywords_begin == '':
            keywords_begin = []
        else:
            keywords_begin = args.keywords_begin.split(",")

        type_seq = args.type_seq
        ignorcase = args.ignorecase
        geo = args.geo
        geo_file = args.geo_file

        species = args.species
        encode_remove = args.encode_remove
        roadmap_remove = args.roadmap_remove

        cwd = args.MetaData
        process = args.process

        if cwd is None:
            cwd = settings['MetaData']

        if cwd is None or cwd == "None":
            cwd = None
            encode_remove = True
            roadmap_remove = True

        SOFTQuickParser(output_prefix, output_path, keywords, keywords_begin, type_seq=type_seq, ignorecase=ignorcase,
                        geo=geo, geofile=geo_file, output_type=species, encode_remove=encode_remove,
                        roadmap_remove=roadmap_remove, encode_pkl=encode_pkl, roadmap_pkl=roadmap_pkl,
                        GSMGSE_pkl=GSMGSE_pkl, cwd=cwd, process=process)
        return

    return 1


def GCF_match():
    '''
    this function provide an entrance to match function

    '''
    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\npython GCF.py match [optional arguments] <first features> <second features>\n\nfor more help, please try: python GCF.py match -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython GCF.py match [optional arguments] <first features> <second features>\n\n",
                                     description='',epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'match' to match samples with different key words within same study")

    parser.add_argument('first_feature_key_words', default=None,
                        help="list of first feature key words need to used to looking for the first set of NGS sequencing samples, "
                             "different key words need to be separated by ','")
    parser.add_argument('second_feature_key_words', default=None,
                        help="list of second feature key words need to used to looking for the second set of NGS sequencing samples, "
                             "different key words need to be separated by ','")

    parser.add_argument('output_prefix1', default=None, help="specify the first output file name prefix.")
    parser.add_argument('output_prefix2', default=None, help="specify the second output file name prefix.")
    parser.add_argument('output_path', default=None, help="specify the output file location.")

    ## optional parameters
    parser.add_argument('-bf', dest='keywords_begin1', default='', metavar='',
                        help="list of key words in first features need to be used to occur in the beginning of a word, "
                             "different key words need to be separated by ','")
    parser.add_argument('-bs', dest='keywords_begin2', default='', metavar='',
                        help="list of key words in second features need to be used to occur in the beginning of a word, "
                             "different key words need to be separated by ','")

    parser.add_argument('-tf', '--typef', dest='type_seq1', default='chip-seq', metavar='',
                        help="first type of sequencing specified in the search. Default is 'chip-seq")
    parser.add_argument('-ts', '--types', dest='type_seq2', default='chip-seq', metavar='',
                        help="second type of sequencing specified in the search. Default is 'chip-seq")

    parser.add_argument('-cf', '--ignorecasef', dest='ignorecase1', default=1, type=int, metavar='',
                        help="specify whether case (A vs a) need to be ignored in the first search. Default is 1 which means case will be ignored in the search. Set to 0 if don't want to ignore the case. ")
    parser.add_argument('-cs', '--ignorecases', dest='ignorecase2', default=1, type=int, metavar='',
                        help="specify whether case (A vs a) need to be ignored in the second search. Default is 1 which means case will be ignored in the search. Set to 0 if don't want to ignore the case. ")

    parser.add_argument('--geo1', dest='geo1', default=0, type=int, metavar='',
                        help="specify whether search will be limited to GEO website search result for the first feature. Default is 0. Set to 1 if want to perform all GEO data search.")
    parser.add_argument('--geo2', dest='geo2', default=0, type=int, metavar='',
                        help="specify whether search will be limited to GEO website search result for the second feature. Default is 0. Set to 1 if want to perform all GEO data search.")


    parser.add_argument('--candidateslistf', dest='geo_file1', default=None, metavar='',
                        help="specify the first file path of GSM ID list if '--hasCandidates' set to 1")
    parser.add_argument('--candidateslists', dest='geo_file2', default=None, metavar='',
                        help="specify the second file path of GSM ID list if '--hasCandidates' set to 1")


    parser.add_argument('-s','--species', dest='species', default='Homo sapiens', metavar='',
                        help="specify the samples' species. Default is Homo sapiens. Please use the species official name. For example, human is Homo sapiens."
                             "If the species name contains space, surround the name with double quotes, for example \"Homo sapiens\"")
    parser.add_argument('-e', '--encode', dest='encode_remove', default=1, type=int, metavar='',
                        help="specify whether need to remove Encode data. Default is 1. Set to 0 to keep Encode data from search.")
    parser.add_argument('-r', '--roadmap', dest='roadmap_remove', default=1, type=int, metavar='',
                        help="specify whether need to remove Roadmap data. Default is 1. Set to 0 to keep Roadmap data from search.")
    parser.add_argument('-m', '--metadata', dest='MetaData', default=None, metavar='',
                        help="specify the GSMs metadata files path")
    parser.add_argument('-p', '--process', dest='process', default=20, type=int, metavar='',
                        help="specify the number of parallel search processes want to use.")

    args = None

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print "\GEO chip seq finder\n"
        parser.print_help()
        print "\n"
        return 0
    elif len(sys.argv) >= 3:
        try:
            args = parser.parse_args()
        except:
            print "\nfor more help, please try: python GCF.py match -h\n"
            return 1

    if args is not None:
        if os.path.exists(args.output_path):
            pass
        else:
            print "Output path is not exist"
            return 1

        settings = get_settings()
        encode_pkl = settings['Encode']
        roadmap_pkl = settings['Roadmap']
        GSMGSE_pkl = settings['GSMGSE_pkl_path']

        keywords1 = args.first_feature_key_words.split(",")
        keywords2 = args.second_feature_key_words.split(",")

        output_prefix1 = args.output_prefix1
        output_prefix2 = args.output_prefix2

        output_path = args.output_path

        if args.keywords_begin1 == '':
            keywords_begin1 = []
        else:
            keywords_begin1 = args.keywords_begin1.split(",")
        if args.keywords_begin2 == '':
            keywords_begin2 = []
        else:
            keywords_begin2 = args.keywords_begin2.split(",")

        type_seq1 = args.type_seq1
        type_seq2 = args.type_seq2

        ignorcase1 = args.ignorecase1
        ignorcase2 = args.ignorecase2

        geo1 = args.geo1
        geo2 = args.geo2

        geo_file1 = args.geo_file1
        geo_file2 = args.geo_file2

        species = args.species

        encode_remove = args.encode_remove
        roadmap_remove = args.roadmap_remove

        cwd = args.MetaData
        process = args.process

        if cwd is None:
            cwd = settings['MetaData']

        if cwd is None or cwd == "None":
            cwd = None
            encode_remove = True
            roadmap_remove = True


        Related_Sample_Search(output_prefix1, output_prefix2, output_path, keywords1, keywords_begin1, keywords2,
                              keywords_begin2,
                              first_type_seq=type_seq1, second_type_seq=type_seq2,
                              first_ignorecase=ignorcase1, second_ignorecase=ignorcase2,
                              first_geo=geo1, first_geofile=geo_file1, second_geo=geo2, second_geofile=geo_file2,
                              output_type=species, encode_remove=encode_remove, roadmap_remove=roadmap_remove,
                              encode_pkl=encode_pkl, roadmap_pkl=roadmap_pkl, GSMGSE_pkl=GSMGSE_pkl, cwd=cwd,
                              process=process)
        return

    return 1


def GCF_query():
    '''
    this function provide an entrance to query function
    '''
    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\npython GCF.py query [optional arguments] <ID list>\n\nfor more help, please try: python GCF.py query -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython GCF.py query [optional arguments] <ID list>\n\n",
                                     description='',epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'query' to looking for samples' NGS sequencing information")

    parser.add_argument('GEO_IDs', default=None,
                        help="list of IDs need to used to looking for the NGS sequencing information, "
                             "it could be a list of IDs separated by ',', or a file containing a list of IDs,"
                             " accepted IDs: GSM, GSE, SRR, SRP, SRX, SAMN, SRP, ",
                        )
    parser.add_argument('output_path', default=None, help="specify the output file name and path.")

    args = None

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print "\GEO chip seq finder\n"
        parser.print_help()
        print "\n"
        return 0
    elif len(sys.argv) >= 3:
        try:
            args = parser.parse_args()
        except:
            print "\nfor more help, please try: python GCF.py query -h\n"
            return 1

    if args is not None:
        if args.output_path.rfind("/") == -1 or os.path.exists(args.output_path[:args.output_path.rfind("/")]):
            pass
        else:
            print "Output path is not exist"

        GEO_ids = args.GEO_IDs
        if os.path.exists(GEO_ids) and os.path.isfile(GEO_ids):
            list_names_obj = open(GEO_ids, "r")
            list_names = []
            for line in list_names_obj.readlines():
                line = line.strip().split(",")
                for l in line:
                    list_names.append(l)
            list_names_obj.close()
            id_list = list(set(list_names))
        else:
            id_list = GEO_ids.split(",")

        settings = get_settings()
        GSMGSE_pkl = settings['GSMGSE_pkl_path']
        GSM_SRR_pkl = settings['GSMtoSRRpkl']

        GEO_query(id_list, args.output_path, GSMGSE_pkl, GSM_SRR_pkl)
        return

    return 1




if len(sys.argv) > 1:
    if sys.argv[1] == "search":
        GCF_search()
    elif sys.argv[1] == "match":
        GCF_match()
    elif sys.argv[1] == "query":
        GCF_query()
    elif sys.argv[1] == "update":
        update()
    elif sys.argv[1] == "setup":
        setup()
    else:
        Help()
else:
    print "\nGEO Chip Finder"
    print "A list of functions for GEO Chip Seq Sample Finder, please try:\npython GCF.py -h"
    print "\nFuctions:"
    print "\tsearch:\n\tsearch chip-seq samples and corresponding input from GEO based on the key words."
    print "\tmatch:\n\tmatch different types of samples from GEO based, for example, looking for the corresponding H3K27me3 samples for each H3K4me3 samples."
    print "\tquery:\n\tget SRR sequencing information by several common identifiers from GEO such as GSE, GSM, SRR, SRX, SRP, etc"
    print ""
