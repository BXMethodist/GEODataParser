import os

output = "./123.txt"

if os.path.exists(output[:output.rfind("/")]):
    pass
else:
    print "Output path is not exist"