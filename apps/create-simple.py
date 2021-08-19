import sys

mat_dir = "/home/shabsheikhha/Documents/BSC/matrices/simple"
lines_num, lines_len = int(sys.argv[1]), int(sys.argv[2])
mat_name = mat_dir + "/simple-" + str(lines_num) + "-" + str(lines_len) + ".mtx"

mat_file = open(mat_name, 'w')

for i in range(lines_num):
    start = i * lines_len
    for j in range(lines_len - 1):
        src = start + j 
        dst = start + (j + 1) % lines_len
        mat_file.write(str(src) + " " + str(dst) + "\n")

mat_file.close()

