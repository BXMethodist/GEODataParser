import os


def get_settings():
    settings = {}
    cwd = os.path.dirname(os.path.realpath(__file__))
    settings_file = open(cwd+'/GCF_settings', "r")
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

    cwd = os.path.dirname(os.path.realpath(__file__))

    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/AIMS/pkl/GSMGSE_map.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/AIMS/pkl/ENCODE_gse.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/AIMS/pkl/Roadmap_gse.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/AIMS/pkl/GGR_gse.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/AIMS/pkl/GSMSRR_map.pkl")
    os.system("wget -P "+cwd+"/pkl/ http://cigwiki.houstonmethodist.org/trackhub/boxia/AIMS/pkl/geoMetaData.db")

    settings = open(cwd + "/GCF_settings", "w")
    settings.write("GSMGSE_pkl_path"+"\t"+cwd+"/pkl/GSMGSE_map.pkl"+"\n")
    settings.write("Encode" + "\t" +cwd+"/pkl/ENCODE_gse.pkl" + "\n")
    settings.write("Roadmap" + "\t" +cwd+"/pkl/Roadmap_gse.pkl" + "\n")
    settings.write("GGR" + "\t" + cwd + "/pkl/GGR_gse.pkl" + "\n")
    settings.write("GSMtoSRRpkl" + "\t" +cwd+"/pkl/GSMSRR_map.pkl" + "\n")
    settings.write("email" + "\t" + "bxia@houstonmethodist.org" + "\n")
    settings.write("MetaData" + "\t" + cwd+"/pkl/geoMetaData.db")
    settings.close()
