from difflib import SequenceMatcher
import re
import numpy as np


a = {}

b = {u'1':'2', u'23':'4'}

b = str(b)[1:-1]

b = b.replace("u'", "")

print b
