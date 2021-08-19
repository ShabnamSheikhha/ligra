import sys

in_name, out_name = sys.argv[1], sys.argv[2]
in_file, out_file = open(in_name, 'r'), open(out_name, 'w')

firstline = True
i = 0
while True:
    line = in_file.readline()
    # reached end of file
    if not line: 
        break
    # skip the comments
    if line[0] == "%":
        continue
    # ignore the first line (contains number of nodes and edges)
    if firstline:
        firstline = False
        continue
    src, dst = map(int, line.split())
    #src -= 1
    #dst -= 1
    out_file.write(str(src) + " " + str(dst) + "\n")
in_file.close()
out_file.close()

