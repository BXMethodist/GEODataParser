import re

a = "CGH08AAA"

b = re.sub("\d", "", a)
print a, b
