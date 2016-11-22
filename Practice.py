
output = open("wgetFiles11-22.txt", "w")
with open("HumanWithH3K4me3Links11-21.txt", "r") as file:
    for line in file.readlines():
        line = line.strip()
        output.write("wget "+ line +';'+'\n')

output.close()
