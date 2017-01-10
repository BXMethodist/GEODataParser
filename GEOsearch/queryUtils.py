import csv
from pickleUtils import *


##/home/tmhbxx3/scratch/XMLhttp/pickles/GSMSRR_map.pkl
###/home/tmhbxx3/scratch/XMLhttp/pickles/SRR_map.pkl

def GSMSRR_query(input_list, output_name, input_type, output_type,
                 GSMtoSRRpkl=None,
                 SRRtoGSMpkl=None):
    if input_type == "GSM" and output_type == "SRR":
        return get_SRR_by_GSM(input_list, output_name, GSMtoSRRpkl, SRRtoGSMpkl)
    elif input_type == "GSM" and output_type == "SRR":
        return get_GSM_by_SRR(input_list, output_name, GSMtoSRRpkl)
    else:
        print "Only information query between GSM and SRR, GSM and GSE is available"

def get_GSM_by_SRR(input_list, output_name, SRRtoGSMpkl):
    if SRRtoGSMpkl is None:
        print "Please provide the SRRtoGSM_map.pkl's location!"
        return

    return


def get_SRR_by_GSM(input_list, output_name, GSMtoSRRpkl, SRRtoGSMpkl):
    '''
    :param input_list: a text file contains a set of input_type files id, each line contains only one id.
    :param output_name: output file path and name.
    :param GSMtoSRRpkl: Location of your GSMSRR_map.pkl
    :param SRRtoGSMpkl: location of your SRR_map.pkl
    :return: a table contains SRRid and GSMid
    '''
    if GSMtoSRRpkl is None or SRRtoGSMpkl is None:
        print "Please provide the SRRtoGSM_map.pkl's location!"
        print "Please provide the GSMtoSRR.pkl's location!"
        return

    GSMtoSRR_map = load_obj(GSMtoSRRpkl)

    SRRtoGSM_map = load_obj(SRRtoGSMpkl)

    result_srrs = set()

    input_file = open(input_list, "r")
    input_ids = [x.strip() for x in input_file.readlines()]
    input_file.close()

    for gsmid in input_ids:
        if gsmid in GSMtoSRR_map:
            result_srrs = result_srrs.union(GSMtoSRR_map[gsmid])

    output_file = open(output_name, "w")
    writer = csv.writer(output_file)
    writer.writerow(["SRR_ID", "GSM_ID", "download_path", "single or paired", "avg_length"])
    for srrid in result_srrs:
        if srrid in SRRtoGSM_map:
            SRR_obj = SRRtoGSM_map[srrid]
            writer.writerow([SRR_obj.SRRid, SRR_obj.gsm_id, SRR_obj.download_path, SRR_obj.layout, SRR_obj.avg_length])
    output_file.close()
    return


if __name__ == "__main__":
    GSMSRR_query("dy_gsm_input.txt", "dy_srr_input.csv", "GSM", "SRR",
                 GSMtoSRRpkl="/home/tmhbxx3/scratch/XMLhttp/pickles/GSMSRR_map.pkl",
                 SRRtoGSMpkl="/home/tmhbxx3/scratch/XMLhttp/pickles/SRR_map.pkl")





