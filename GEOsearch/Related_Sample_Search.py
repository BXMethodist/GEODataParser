from GEOsearch import SOFTQuickParser
from collections import defaultdict
from input_search_utils import keyword, Similarity, Character_Similarity


def Related_Sample_Search(output_surfix, first_features, first_features_begin, second_features, second_begin_features,
                    first_type_seq="chip-seq", second_type_seq="chip-seq", cwd=None, first_ignorecase=True, second_ignorecase=True, geo=False, first_geofile=None, second_geofile=None, output_type="Homo sapiens",
                    encode_remove=False, roadmap_remove=False):
    first_samples = SOFTQuickParser(output_surfix, first_features, first_features_begin, type_seq=first_type_seq,
                                    cwd=cwd, ignorecase=first_ignorecase, geo=geo, geofile=first_geofile,
                                    output_type=output_type, encode_remove=encode_remove, roadmap_remove=roadmap_remove)
    second_samples = SOFTQuickParser(output_surfix, second_features, second_begin_features, type_seq=second_type_seq,
                                     cwd=cwd, ignorecase=second_ignorecase, geo=geo, geofile=second_geofile,
                                     output_type=output_type, encode_remove=encode_remove, roadmap_remove=roadmap_remove)

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
                related_feature_key_word = keyword(related_sample.title, second_features, second_begin_features, second_ignorecase)

                if sample.cellType == related_sample.cellType:
                    score = Similarity(sample.title, feature_key_word, related_sample.title, related_feature_key_word)

                if score > best_score:
                    best_score = score
                    best_id = set()
                    best_id.add(related_sample.id)
                elif score == best_score:
                    best_id.add(related_sample.id)

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
            targetGSMs = targetGSMs.union(second_groupbygse[gse])
        if sample.title_found:
            feature_key_word = keyword(sample.title, first_features, first_features_begin, first_ignorecase)

            best_score = float("-inf")
            best_id = set()

            for related_id in targetGSMs:
                score = None
                related_sample = first_samples[related_id]
                related_feature_key_word = keyword(related_sample.title, second_features, second_begin_features, second_ignorecase)

                if sample.cellType == related_sample.cellType:
                    score = Similarity(sample.title, feature_key_word, related_sample.title, related_feature_key_word)

                if score > best_score:
                    best_score = score
                    best_id = set()
                    best_id.add(related_sample.id)
                elif score == best_score:
                    best_id.add(related_sample.id)

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

    ### cross validation
    for key, values in pairs.items():
        for value in values:
            if key not in pairs[value]:
                print key, value, "cross validation failed!!!!"

    #### save file
    import csv
    output = open(output_surfix, "w")
    writer = csv.writer(output)
    writer.writerow(["sample1_id", "sample1_descriptions", "sample1_series_id", "sample2_id", "sample2_descriptions", "sample2_series_id"])
    for key, values in pairs.items():
        for value in values:
            writer.writerow([key, first_samples[key].features, first_samples[key].series,
                             value, second_samples[key].features, second_samples[key].series])
    output.close()
    return pairs

if __name__ == "__main__":
    Related_Sample_Search("H3K4me3_vs_H3K27me1", ["h3k4me3", "k4me3", "k4m3", "h3k4m3"],[],
                          ["h3k27me1", "k27me1", "k27m1","h3k27m1"], [],
                          cwd="/home/tmhbxx3/scratch/XMLhttp/QuickXMLs")
