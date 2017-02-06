from difflib import SequenceMatcher
import re


def Similarity(title1, keyword1, title2, keyword2):
    score = SequenceMatcher(None, title1, title2).ratio()

    title1 = title1.lower().replace(keyword1.lower(), "").replace("chip-seq", "")
    title2 = title2.lower().replace(keyword2.lower(), "").replace("chip-seq", "")

    title1 = re.sub(r'rep[0-9]*', '', title1)
    title2 = re.sub(r'rep[0-9]*', '', title2)

    print title1, title2

    score_replace = SequenceMatcher(None, title1, title2).ratio()

    print score_replace

    return max(score, score_replace)


a = 'MCF10A_H3K4ME3_REP1'
b = 'MCF10A_H3K4ME3_REP2'
c = 'MCF10A_input_REP1'
d = 'MCF10A_input_REP2'

print Similarity(a, 'H3K4me3', c, 'input')
print Similarity(a, 'H3K4me3', d, 'input')

print Similarity(a, 'H3K4me3', c, 'input') == Similarity(a, 'H3K4me3', d, 'input')
