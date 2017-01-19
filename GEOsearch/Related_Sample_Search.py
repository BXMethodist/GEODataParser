from GEOsearch import SOFTQuickParser
from collections import defaultdict
from input_search_utils import keyword, Similarity, Character_Similarity


def Related_Sample_Search(output_surfix1, output_surfix2,  first_features, first_features_begin, second_features, second_begin_features,
                    first_type_seq="chip-seq", second_type_seq="chip-seq", first_ignorecase=True, second_ignorecase=True,
                    first_geo=False, first_geofile=None, second_geo=False, second_geofile=None, output_type="Homo sapiens",
                    encode_remove=False, roadmap_remove=False, encode_pkl=None, roadmap_pkl=None, GSMGSE_pkl=None):

    first_samples = SOFTQuickParser(output_surfix1, first_features, first_features_begin, type_seq=first_type_seq,
                                    ignorecase=first_ignorecase, geo=first_geo, geofile=first_geofile,
                                    output_type=output_type, encode_remove=encode_remove, roadmap_remove=roadmap_remove,
                                    encode_pkl=encode_pkl, roadmap_pkl=roadmap_pkl, GSMGSE_pkl=GSMGSE_pkl)

    second_samples = SOFTQuickParser(output_surfix2, second_features, second_begin_features, type_seq=second_type_seq,
                                     ignorecase=second_ignorecase, geo=second_geo, geofile=second_geofile,
                                     output_type=output_type, encode_remove=encode_remove, roadmap_remove=roadmap_remove,
                                     encode_pkl=encode_pkl, roadmap_pkl=roadmap_pkl, GSMGSE_pkl=GSMGSE_pkl)

    for id in first_samples.keys():
        if id in second_samples:
            print id, " in both results"



    ###TODO: find common gse and decide which samples are a pair
    ## STEP1 groupbygse
    first_groupbygse = defaultdict(set)

    for sample in first_samples.values():
        for gse in sample.series:
            first_groupbygse[gse].add(sample.id)

    second_groupbygse = defaultdict(set)

    for sample in second_samples.values():
        for gse in sample.series:
            second_groupbygse[gse].add(sample.id)

    pairs = defaultdict(set)

    for sample in first_samples.values():
        targetGSMs = set()
        for gse in sample.series:
            targetGSMs = targetGSMs.union(second_groupbygse[gse])
        if sample.title_found:
            feature_key_word = keyword(sample.title, first_features, first_features_begin, first_ignorecase)

            best_score = float("-inf")
            best_id = set()

            for related_id in targetGSMs:
                score = None
                related_sample = second_samples[related_id]
                if related_sample.title_found:
                    related_feature_key_word = keyword(related_sample.title, second_features, second_begin_features, second_ignorecase)

                    if sample.cellType == related_sample.cellType:
                        score = Similarity(sample.title, feature_key_word, related_sample.title, related_feature_key_word)

                    if score > best_score:
                        best_score = score
                        best_id = set()
                        best_id.add(related_sample.id)
                    elif score == best_score:
                        best_id.add(related_sample.id)
                else:
                    continue

            if best_id:
                pairs[sample.id] = pairs[sample.id].union(best_id)

        else:
            # TODO none-title sample
            best_score = float("-inf")
            best_id = set()
            for related_id in targetGSMs:
                score = None
                related_sample = second_samples[related_id]

                if sample.cellType == related_sample.cellType:
                    score = Character_Similarity(sample, related_sample)

                if score > best_score:
                    best_score = score
                    best_id = set()
                    best_id.add(related_sample.id)
                elif score == best_score:
                    best_id.add(related_sample.id)
            if best_id:
                pairs[sample.id] = pairs[sample.id].union(best_id)

    for sample in second_samples.values():
        targetGSMs = set()
        for gse in sample.series:
            targetGSMs = targetGSMs.union(first_groupbygse[gse])
        if sample.title_found:
            feature_key_word = keyword(sample.title, second_features, second_begin_features, second_ignorecase)

            best_score = float("-inf")
            best_id = set()

            for related_id in targetGSMs:
                score = None
                related_sample = first_samples[related_id]
                if related_sample.title_found:
                    related_feature_key_word = keyword(related_sample.title, first_features, first_features_begin, first_ignorecase)

                    if sample.cellType == related_sample.cellType:
                        score = Similarity(sample.title, feature_key_word, related_sample.title, related_feature_key_word)

                    if score > best_score:
                        best_score = score
                        best_id = set()
                        best_id.add(related_sample.id)
                    elif score == best_score:
                        best_id.add(related_sample.id)
                else:
                    continue

            if best_id:
                pairs[sample.id] = pairs[sample.id].union(best_id)

        else:
            # TODO none-title sample
            best_score = float("-inf")
            best_id = set()
            for related_id in targetGSMs:
                score = None
                related_sample = first_samples[related_id]

                if sample.cellType == related_sample.cellType:
                    score = Character_Similarity(sample, related_sample)

                if score > best_score:
                    best_score = score
                    best_id = set()
                    best_id.add(related_sample.id)
                elif score == best_score:
                    best_id.add(related_sample.id)
            if best_id:
                pairs[sample.id] = pairs[sample.id].union(best_id)


    wrong_pairs = set()
    new_pairs = set()
    ### cross validation
    for key, values in pairs.items():
        for value in values:
            if key not in pairs[value]:
                wrong_pairs.add((key, value))
            elif (key, value) not in new_pairs and (value, key) not in new_pairs:
                new_pairs.add((key, value))

    #### save file
    import csv
    output = open(output_surfix1+output_surfix2+".txt", "w")
    writer = csv.writer(output)
    writer.writerow(["sample1_id", "sample1_title", "sample1_descriptions", "sample1_series_id", "sample2_id", "sample2_title", "sample2_descriptions", "sample2_series_id"])
    for pair in new_pairs:
        key, value = pair
        if key in first_samples:
            try:
                writer.writerow([key, first_samples[key].title, first_samples[key].features, first_samples[key].series,
                                 value, second_samples[value].title, second_samples[value].features, second_samples[value].series])
            except:
                writer.writerow([value, first_samples[value].title, first_samples[value].features, first_samples[value].series,
                                 key, second_samples[key].title, second_samples[key].features, second_samples[key].series,])
        elif key in second_samples:
            try:
                writer.writerow([value, first_samples[value].title, first_samples[value].features, first_samples[value].series,
                                 key, second_samples[key].title, second_samples[key].features, second_samples[key].series,
                     ])
            except:
                writer.writerow([key, first_samples[key].title, first_samples[key].features, first_samples[key].series,
                     value, second_samples[value].title, second_samples[value].features, second_samples[value].series])
        else:
            print key
    output.close()
    return pairs

if __name__ == "__main__":
    Related_Sample_Search("H3K4me3", "H3K27me3", ["h3k4me3", "k4me3", "k4m3", "h3k4m3"],[],
                          ["h3k27me3", "k27me3", "k27m3","h3k27m3"], [],
                          first_geo=False, first_geofile="uniqueGSM_GEOsearch.txt",
                          second_geo=False, second_geofile="unique_H3K27me3.txt", encode_remove=True, roadmap_remove=True)
