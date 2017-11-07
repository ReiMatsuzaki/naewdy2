import pandas as pd
from naewdy2.math import ijkv2ten

df = pd.DataFrame({"i":[1,1,3], "j":[1,2,3], "k":[2,1,2], "val":[1.1,1.2,1.3]})
ten = ijkv2ten(df)
print ten








