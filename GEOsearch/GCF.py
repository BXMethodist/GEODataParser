import argparse, sys, os
from queryUtils import GEO_query


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
    pass


def GCF_match():
    '''
    this function provide an entrance to match function

    '''
    pass


def GCF_query():
    '''
    this function provide an entrance to query function
    '''
    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\npython GCF.py query [optional arguments] <ID list>\n\nfor more help, please try: python danpos.py query -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython GCF.py  <command> [optional arguments] <ID list>\n\n",
                                     description='',epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'query' to looking for samples' NGS sequencing information")

    parser.add_argument('GEO_IDs', default=None,
                        help="list of IDs need to used to looking for the NGS sequencing information, "
                             "it could be a list of IDs separated by ',', or a file containing a list of IDs")
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
        GEO_ids = args.GEO_IDs
        if os.path.exists(GEO_ids) and os.path.isfile(GEO_ids):
            list_names_obj = open(GEO_ids, "r")
            list_names = [x.strip() for x in list_names_obj.readlines()]
            list_names_obj.close()
            id_list = list_names
        else:
            id_list = GEO_ids.split(",")

        GEO_query(id_list, args.output_path)
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
