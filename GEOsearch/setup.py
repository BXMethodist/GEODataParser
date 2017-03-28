"""
Copyright (c) <2017> <Dr. Kaifu Chen lab, Research Institute, Houston Methodist Hospital >

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

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
