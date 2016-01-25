f = open('test.prm','r')
f2 = open('test.ppm','w')
o = 0
columns = 0
for line in f:
    l = 0
    line2 = line.split(" ")
    o = o + 1
    filelength = (len(line2)-2)/3 
    if o > 1:
        if o == 2:
            columns = int(line)
            f2.write("P3\n" + str(columns) + " " + str(columns) + "\n255\n")
        else:
            for p in line2:
                l = l + 1
                if l/3 == columns*columns/2:
                    f2.write(line2[o] + " " + line2[o + 1] + " " + line2[o + 2] + " ")
f.close()
f2.close()
