from difflib import SequenceMatcher
import re

def isInput(title, feature_key_word):
    non_capital_keywords = ['input','wce']

    if feature_key_word.find("H3K") != -1:
        capital_keywords = ['IgG']
    else:
        capital_keywords = ['IgG', '_H3_', " H3"]

    for c in non_capital_keywords:
        if title.lower().find(c) != -1:
            return True, c

    for n in capital_keywords:
        if n ==' H3' and title.endswith(n):
            return True, 'H3'
        elif title.find(n) != -1:
            if n == "_H3_":
                return True, 'H3'
            return True, n
    return False, ""


def has_antibody(sample, keyword):
    for v in sample.antibody.values():
        if v.find(keyword) != -1:
            return True
    return False

print isInput('Exp015-2_H3K4me3_500c_Ad_BC01_300u', 'H3K4me3')
