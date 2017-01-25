import numpy as np, pandas as pd

a = np.arange(100)

df = pd.DataFrame(a)

df.to_csv("./123/123.csv")