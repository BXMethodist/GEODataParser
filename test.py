import pandas as pd


a = [[1,2,3],[4,5,6]]
b = [[1,2],[3,4]]

dfa = pd.DataFrame(a, columns=['A','B','C'])
dfb = pd.DataFrame(b, columns=['A', 'B'])

print dfa.append(dfb)