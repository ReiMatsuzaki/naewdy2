from naewdy2.io import *

with open("tmp0.dat", "w") as f:
    f.write("""aa
$VEC
11A1
1B11
C111
$END
other comment
hhhhhhhhhhhhhhhhhhh
$VEC
222A
2B22
$END
cococooc
jijiji
jkajfd
$VEC
A333
33B3
$END
finish
""")

with open("tmp1.dat", "w") as f:
    f.write("""hello
vec:
""")
extract_lines("$VEC", "$END", "tmp0.dat", "tmp1.dat", "a", 1)
extract_lines("$VEC", "$END", "tmp0.dat", "tmp1.dat", "a", 2)

with open("tmp2.dat", "w") as f:
    f.write("""Hello
__x__
Hello?
__y__
""")
    
rep_copy([("__x__", "strx"), ("__y__", "stry")], "tmp2.dat", "tmp3.dat")
             
