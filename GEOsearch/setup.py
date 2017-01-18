import os


def setup(command):
    os.system("wget -P ./pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/GSMGSE_map.pkl")
    os.system("wget -P ./pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/ENCODE_gse.pkl")
    os.system("wget -P ./pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/Roadmap_gse.pkl")
    os.system("wget -P ./pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/GSMSRR_map.pkl")

    settings = open("GCF_settings.txt", "w")
    settings.write("GSMGSE_pkl_path"+"\t"+"./pkl/GSMGSE_map.pkl"+"\n")
    settings.write("Encode" + "\t" + "./pkl/ENCODE_gse.pkl" + "\n")
    settings.write("Roadmap" + "\t" + "./pkl/Roadmap_gse.pkl" + "\n")
    settings.write("GSMtoSRRpkl" + "\t" + "./pkl/GSMSRR_map.pkl" + "\n")

    if command == "query":
        settings.close()
        return 0
    elif command == "all":
        settings.write("MetaData" + "\t" + "./MetaData/GSM" + "\n")
        settings.close()
        os.system("wget -P ./MetaData/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/MetaData/GSE.tar.gz")
        os.system("mkdir ./MetaData/GSE")
        os.system("tar -xvzf ./MetaData/GSE.tar.gz -C ./MetaData/GSE")
        os.system("wget -P ./MetaData/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/MetaData/GSM.tar.gz")
        os.system("mkdir ./MetaData/GSM")
        os.system("tar -xvzf ./MetaData/GSM.tar.gz -C ./MetaData/GSM")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython setup.py <command>\n\n",
                                     description='', epilog="Chen lab, Houston Methodist")

    parser.add_argument('command', default=None, help="set as 'query' if only need to use NGS sample information query function,"
                                                      "set as 'all' if need search and match NGS sample with key words functions")