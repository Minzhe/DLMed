import sys

path = sys.argv[1]
res = []
with open(path, 'r') as f:
    flag = -1
    for line in f:
        if line.strip().startswith('***'):
            flag *= -1
        if flag == 1:
            res.append(line.strip())
        
for line in res:
    print(line)