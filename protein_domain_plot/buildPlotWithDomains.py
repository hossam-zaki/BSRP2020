import re
import subprocess

result = subprocess.run(
    ["Rscript", "protein_domain_plot/getDomains.R", "--start", "500", "--width", "30", "--transcript", "ENST00000487861"], stdout=subprocess.PIPE, text=True)

str = result.stdout.split('\n')

start = []
for sub in str:
    start.append(re.match(".*ENSP00000419881\s+(\d+)\s+(\d+).*", sub))

for i in start:
    if i:
        print("yee")
        print(i[1])
        print("yee2")
        print(i[2])
