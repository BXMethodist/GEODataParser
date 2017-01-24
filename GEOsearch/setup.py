import os


def get_settings():
    settings = {}
    settings_file = open('GCF_settings', "r")
    for line in settings_file.readlines():
        info = line.split()
        settings[info[0]] = info[1].strip()
    settings_file.close()
    return settings


def setup():
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython setup.py <command>\n\n",
                                     description='', epilog="Chen lab, Houston Methodist")

    parser.add_argument('email', default=None, help="email is required to connect to NCBI ftp site")

    args = parser.parse_args()

    email = args.email

    cwd = os.getcwd()

    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/GSMGSE_map.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/ENCODE_gse.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/Roadmap_gse.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/GSMSRR_map.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/GCF/pkl/geoMetaData.db")

    settings = open("GCF_settings.txt", "w")
    settings.write("GSMGSE_pkl_path"+"\t"+cwd+"/pkl/GSMGSE_map.pkl"+"\n")
    settings.write("Encode" + "\t" +cwd+"/pkl/ENCODE_gse.pkl" + "\n")
    settings.write("Roadmap" + "\t" +cwd+"/pkl/Roadmap_gse.pkl" + "\n")
    settings.write("GSMtoSRRpkl" + "\t" +cwd+"/pkl/GSMSRR_map.pkl" + "\n")
    settings.write("email" + "\t" + email + "\n")
    settings.write("MetaData" + "\t" + cwd+"/pkl/geoMetaData.db")
    settings.close()
