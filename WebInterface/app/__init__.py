from flask import Flask, render_template, request
from queryUtils import GEO_query
from search import SOFTQuickParser
from Related_Sample_Search import Related_Sample_Search
from setup import get_settings, setup
from update import update
import multiprocessing

app = Flask(__name__)

@app.route("/")
def main():
    return render_template('index.html')

@app.route('/search')
def showSearch():
    return render_template('search.html')

@app.route('/match')
def showMatch():
    return render_template('match.html')

@app.route('/query')
def showQuery():
    return render_template('query.html')

@app.route('/search',methods=['POST'])
def search():
    _searchterms = request.form['searchterms']
    _Species = request.form['Species']
    _inputEmail = request.form['inputEmail']

    settings = get_settings()
    encode_pkl = settings['Encode']
    roadmap_pkl = settings['Roadmap']
    GGRmap_pkl = settings['GGR']
    GSMGSE_pkl = settings['GSMGSE_pkl_path']

    keywords = _searchterms.split(",")

    output_prefix = keywords[0]

    output_path = '../tmp/'

    keywords_begin = []

    type_seq = 'chip-seq'
    ignorcase = True
    geo = False
    geo_file = None

    species = _Species
    encode_remove = True
    roadmap_remove = True

    cwd = settings['MetaData']
    process = 20

    pool = multiprocessing.Pool(processes=1)
    pool.apply_async(SOFTQuickParser, args=(output_prefix, output_path, output_path,
                                            keywords, keywords_begin, type_seq,
                                            ignorcase, geo, geo_file, species, encode_remove, roadmap_remove,
                                            encode_pkl, roadmap_pkl, GGRmap_pkl,
                                            GSMGSE_pkl, cwd, process, _inputEmail))
    return 'We are processing your request, results will be sent to your email'

@app.route('/match',methods=['POST'])
def match():
    _feature1 = request.form['feature1']
    _feature2 = request.form['feature2']
    _Species = request.form['Species']

    _inputEmail = request.form['inputEmail']

    settings = get_settings()
    encode_pkl = settings['Encode']
    roadmap_pkl = settings['Roadmap']
    GGRmap_pkl = settings['GGR']
    GSMGSE_pkl = settings['GSMGSE_pkl_path']

    keywords1 = _feature1.split(",")
    keywords2 = _feature2.split(",")
    output_prefix1 = keywords1[0]
    output_prefix2 = keywords2[0]
    output_path = '../tmp/'
    type_seq1 = 'chip-seq'
    type_seq2 = 'chip-seq'

    species = _Species if _Species != '' else 'Homo sapiens'
    cwd = settings['MetaData']

    pool = multiprocessing.Pool(processes=1)
    pool.apply_async(Related_Sample_Search, args=(output_prefix1, output_prefix2, output_path,
                          keywords1, [], keywords2, [],
                          type_seq1, type_seq2, True, True, False, None, False, None,
                          species, True, True,
                          encode_pkl, roadmap_pkl, GGRmap_pkl,
                          GSMGSE_pkl, cwd, 20, _inputEmail))
    return 'We are processing your request, results will be sent to your email'

@app.route('/query',methods=['GET', 'POST'])
def query():
    f = request.files['IDlist']
    id_list = []
    for line in f.readlines():
        id_list.append(line.strip())
    f.close()
    output_path = '../tmp/query.txt'
    settings = get_settings()
    GSMGSE_pkl = settings['GSMGSE_pkl_path']
    GSM_SRR_pkl = settings['GSMtoSRRpkl']

    pool = multiprocessing.Pool(processes=1)
    pool.apply_async(GEO_query, args=(id_list, output_path, GSMGSE_pkl, GSM_SRR_pkl))
    return 'We are processing your request, results will be sent to your email'

if __name__ == "__main__":
    app.run()


# from flask import Flask, render_template, request
#
# app = Flask(__name__)
#
#
# @app.route('/')
# def upload():
#     return render_template('upload.html')
#
#
# @app.route('/upload', methods=['GET', 'POST'])
# def upload_file():
#     if request.method == 'POST':
#         f = request.files['file']
#         f.save(f.filename)
#         return 'file uploaded successfully'
#
#
# if __name__ == '__main__':
#     app.run()