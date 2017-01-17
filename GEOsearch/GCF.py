import argparse, sys, os
from queryUtils import GEO_query
from GEOsearch import SOFTQuickParser
from Related_Sample_Search import Related_Sample_Search


def Help():
    print "\nGEO Chip Finder"
    print "A list of functions for GEO Chip Seq Sample Finder, please try:\npython GCF.py -h"
    print "\nFuctions:"
    print "\tsearch:\n\tsearch chip-seq samples and corresponding input from GEO based on the key words."
    print "\tmatch:\n\tmatch different types of samples from GEO based, for example, looking for the corresponding H3K27me3 samples for each H3K4me3 samples."
    print "\tquery:\n\tget SRR sequencing information by several common identifiers from GEO such as GSE, GSM, SRR, SRX, SRP, etc"


def get_settings():
    settings = {}
    settings_file = open('GCF_settings', "r")
    for line in settings_file.readlines():
        info = line.split()
        settings[info[0]] = info[1].strip()
    settings_file.close()
    return settings


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

    parser.add_argument('feature key words', dest='keywords', default=None,
                        help="list of feature key words need to used to looking for the NGS sequencing samples, "
                             "different key words need to be separated by ','")
    parser.add_argument('output_prefix', default=None, help="specify the output file prefix.")

    ## optional parameters
    parser.add_argument('-b', dest='keywords_begin', default='',
                        help="list of feature key words need to used to occur in the beginning of a word, "
                             "different key words need to be separated by ','")
    parser.add_argument('-t', '--type', dest='type_seq', default='chip-seq', help="type of sequencing specified in the search. Default is 'chip-seq")
    parser.add_argument('-c', '--ignorecase', dest='ignorecase', default=1, type=int,
                        help="specify whether case (A vs a) need to be ignored in the search. Default is 1 which means case will be ignored in the search. Set to 0 if don't want to ignore the case. ")
    parser.add_argument('--hasCandidates', dest='geo', default=0,
                        help="specify whether there will be a pre-search GSM ID list provided. Default is 0. Set to 1 if GSM ID list will be provided.")
    parser.add_argument('--candidateslist', dest='geo_file', default=None,
                        help="specify the file path of GSM ID list if '--hasCandidates' set to 1")
    parser.add_argument('-s','--species', dest='species', default='Homo sapiens',
                        help="specify the samples' species. Default is Homo sapiens. Please use the species official name. For example, human is Homo sapiens.")
    parser.add_argument('-e', '--encode', dest='encode_remove', default=0, type=int,
                        help="specify whether need to remove Encode data. Default is 0. Set to 1 to remove Encode data from search.")
    parser.add_argument('-r', '--roadmap', dest='roadmap_remove', default=0, type=int,
                        help="specify whether need to remove Roadmap data. Default is 0. Set to 1 to remove Roadmap data from search.")

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
        settings = get_settings()
        metadata_path = settings['MetaData']
        encode_pkl = settings['Encode']
        roadmap_pkl = settings['Roadmap']
        GSMGSE_pkl = settings['GSMGSE_pkl_path']

        keywords = args.keywords.split(",")
        output_path = args.output_prefix
        keywords_begin = args.keywords_begin.split(",")
        type_seq = args.type_seq
        ignorcase = args.ignorecase
        geo = args.geo
        geo_file = args.geo_file

        if geo and geo_file is None:
            print "Please provide the file location for candidates GSM list."
            return 1

        species = args.species
        encode_remove = args.encode_remove
        roadmap_remove = args.roadmap_remove

        SOFTQuickParser(output_path, keywords, keywords_begin, type_seq=type_seq, cwd=metadata_path, ignorecase=ignorcase,
                        geo=geo, geofile=geo_file, output_type=species, encode_remove=encode_remove,
                        roadmap_remove=roadmap_remove, encode_pkl=encode_pkl, roadmap_pkl=roadmap_pkl,
                        GSMGSE_pkl=GSMGSE_pkl)
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

    parser.add_argument('first feature key words', dest='keywords1', default=None,
                        help="list of first feature key words need to used to looking for the first set of NGS sequencing samples, "
                             "different key words need to be separated by ','")
    parser.add_argument('second feature key words', dest='keywords2', default=None,
                        help="list of second feature key words need to used to looking for the second set of NGS sequencing samples, "
                             "different key words need to be separated by ','")

    parser.add_argument('output_prefix1', default=None, help="specify the first output file name prefix.")
    parser.add_argument('output_prefix2', default=None, help="specify the second output file name prefix.")

    ## optional parameters
    parser.add_argument('-bf', dest='keywords_begin1', default='',
                        help="list of first feature key words need to used to occur in the beginning of a word, "
                             "different key words need to be separated by ','")
    parser.add_argument('-bs', dest='keywords_begin2', default='',
                        help="list of second feature key words need to used to occur in the beginning of a word, "
                             "different key words need to be separated by ','")

    parser.add_argument('-tf', '--typef', dest='type_seq1', default='chip-seq',
                        help="first type of sequencing specified in the search. Default is 'chip-seq")
    parser.add_argument('-ts', '--types', dest='type_seq2', default='chip-seq',
                        help="second type of sequencing specified in the search. Default is 'chip-seq")

    parser.add_argument('-cf', '--ignorecasef', dest='ignorecase1', default=1, type=int,
                        help="specify whether case (A vs a) need to be ignored in the first search. Default is 1 which means case will be ignored in the search. Set to 0 if don't want to ignore the case. ")
    parser.add_argument('-cs', '--ignorecases', dest='ignorecase2', default=1, type=int,
                        help="specify whether case (A vs a) need to be ignored in the second search. Default is 1 which means case will be ignored in the search. Set to 0 if don't want to ignore the case. ")

    parser.add_argument('--hasCandidatesf', dest='geo1', default=0,
                        help="specify whether there will be a pre-search GSM ID list provided for the first search. Default is 0. Set to 1 if GSM ID list will be provided.")
    parser.add_argument('--hasCandidatess', dest='geo2', default=0,
                        help="specify whether there will be a pre-search GSM ID list provided for the second search. Default is 0. Set to 1 if GSM ID list will be provided.")


    parser.add_argument('--candidateslistf', dest='geo_file1', default=None,
                        help="specify the first file path of GSM ID list if '--hasCandidates' set to 1")
    parser.add_argument('--candidateslists', dest='geo_file2', default=None,
                        help="specify the second file path of GSM ID list if '--hasCandidates' set to 1")


    parser.add_argument('-s','--species', dest='species', default='Homo sapiens',
                        help="specify the samples' species. Default is Homo sapiens. Please use the species official name. For example, human is Homo sapiens.")
    parser.add_argument('-e', '--encode', dest='encode_remove', default=0, type=int,
                        help="specify whether need to remove Encode data. Default is 0. Set to 1 to remove Encode data from search.")
    parser.add_argument('-r', '--roadmap', dest='roadmap_remove', default=0, type=int,
                        help="specify whether need to remove Roadmap data. Default is 0. Set to 1 to remove Roadmap data from search.")

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
        settings = get_settings()
        metadata_path = settings['MetaData']
        encode_pkl = settings['Encode']
        roadmap_pkl = settings['Roadmap']
        GSMGSE_pkl = settings['GSMGSE_pkl_path']

        keywords1 = args.keywords1.split(",")
        keywords2 = args.keywords2.split(",")

        output_path1 = args.output_prefix1
        output_path2 = args.output_prefix2

        keywords_begin1 = args.keywords_begin1.split(",")
        keywords_begin2 = args.keywords_begin2.split(",")

        type_seq1 = args.type_seq1
        type_seq2 = args.type_seq2

        ignorcase1 = args.ignorecase1
        ignorcase2 = args.ignorecase2

        geo1 = args.geo1
        geo2 = args.geo2

        geo_file1 = args.geo_file1
        geo_file2 = args.geo_file2

        if geo1 and geo_file1 is None:
            print "Please provide the file location for the first candidates GSM list."
            return 1
        if geo2 and geo_file2 is None:
            print "Please provide the file location for the second candidates GSM list."
            return 1

        species = args.species
        encode_remove = args.encode_remove
        roadmap_remove = args.roadmap_remove

        Related_Sample_Search(output_path1, output_path2, keywords1, keywords_begin1, keywords2,
                              keywords_begin2,
                              first_type_seq=type_seq1, second_type_seq=type_seq2, cwd=metadata_path,
                              first_ignorecase=ignorcase1, second_ignorecase=ignorcase2,
                              first_geo=geo1, first_geofile=geo_file1, second_geo=geo2, second_geofile=geo_file2,
                              output_type=species, encode_remove=encode_remove, roadmap_remove=roadmap_remove,
                              encode_pkl=encode_pkl, roadmap_pkl=roadmap_pkl, GSMGSE_pkl=GSMGSE_pkl)
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
                             "it could be a list of IDs separated by ',', or a file containing a list of IDs")
    parser.add_argument('output_path', default=None, help="specify the output file name and path.")

    parser.add_argument('-g', '--gse', dest='hasGSE', default=0, type=int, help="indicate whether query contains GSE ID, default is 0 (False). If query contains GSE ID, set -g to 1.")

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
        GEO_ids = args.GEO_IDs
        if os.path.exists(GEO_ids) and os.path.isfile(GEO_ids):
            list_names_obj = open(GEO_ids, "r")
            list_names = [x.strip() for x in list_names_obj.readlines()]
            list_names_obj.close()
            id_list = list_names
        else:
            id_list = GEO_ids.split(",")

        if not args.hasGSE:
            GSMGSE_pkl_path = None
        else:
            settings = get_settings()
            GSMGSE_pkl_path = settings['GSMGSE_pkl_path']

        GEO_query(id_list, args.output_path, GSMGSE_pkl_path, args.hasGSE)
        return

    return 1



if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "search":
            GCF_search()
        elif sys.argv[1] == "match":
            GCF_match()
        elif sys.argv[1] == "query":
            GCF_query()
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
