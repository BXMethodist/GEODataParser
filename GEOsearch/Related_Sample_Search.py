from GEOsearch import SOFTQuickParser

if __name__ == '__main__':
    def Related_Sample_Search(output_surfix, first_features, first_features_begin, second_features, second_begin_features,
                        first_type_seq="chip-seq", second_type_seq="chip-seq", cwd=None, first_ignorecase=True, second_ignorecase=True, geo=False, first_geofile=None, second_geofile=None, output_type="Homo sapiens",
                        encode_remove=False, roadmap_remove=False):
        first_samples = SOFTQuickParser(output_surfix, first_features, first_features_begin, type_seq=first_type_seq, cwd=cwd, ignorecase=first_ignorecase, geo=geo, geofile=first_geofile, output_type=output_type)
        second_samples = SOFTQuickParser(output_surfix, second_features, second_begin_features, type_seq=second_type_seq,
                                         cwd=cwd, ignorecase=second_ignorecase, geo=geo, geofile=second_geofile,
                                         output_type=output_type)

        ###TODO: find common gse and decide which samples are a pair

        return
