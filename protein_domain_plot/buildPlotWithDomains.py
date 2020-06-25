import re
import subprocess

result = subprocess.run(
    ["Rscript", "protein_domain_plot/getDomains.R", "--start", "500", "--width", "50", "--transcript", "ENST00000487861"], stdout=subprocess.PIPE, text=True)

str = result.stdout.split('\n')

start = []
bool = False
for sub in str:
    print(bool)
    if bool:
        start.append(re.match(".+\s+(\d+).+", sub)[1])
        bool = False
    if '<integer> <logical>' in sub:
        bool = True
        continue
    match = re.match(".+ENST00000487861\s+(\d+).*", sub)
    if match:
        start.append(match[1])
for i in start:
    if i:
        print(i)
