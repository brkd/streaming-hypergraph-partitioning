import matplotlib.pyplot as plt

X, Y_filler = [], []
with open("N2Pvertex.txt") as n2pfile:
    ctr = 0
    for line in n2pfile:
        if '%' in line:
            ctr += 1
            if ctr != 2:
                continue
            else:
                break
        X.append(int(line.split(',')[0]))
        Y_filler.append(int(line.split(',')[2]))    

Y = []
for i in range(0, len(Y_filler)):
    if i == 0:
        Y.append(Y_filler[i])
    else:
        Y.append(Y[i - 1] + Y_filler[i])

X2, Y_filler2 = [], []
with open("BFvertex.txt") as bffile:
    ctr = 0
    for line in bffile:
        if '%' in line:
            ctr += 1
            if ctr != 2:
                continue
            else:
                break
        X2.append(int(line.split(',')[0]))
        Y_filler2.append(int(line.split(',')[2]))    

Y2 = []
for i in range(0, len(Y_filler2)):
    if i == 0:
        Y2.append(Y_filler2[i])
    else:
        Y2.append(Y2[i - 1] + Y_filler2[i])

X3, Y_filler3 = [], []
with open("BF2vertex.txt") as bf2file:
    ctr = 0
    for line in bf2file:
        if '%' in line:
            ctr += 1
            if ctr != 2:
                continue
            else:
                break
        X3.append(int(line.split(',')[0]))
        Y_filler3.append(int(line.split(',')[2]))    

Y3 = []
for i in range(0, len(Y_filler3)):
    if i == 0:
        Y3.append(Y_filler3[i])
    else:
        Y3.append(Y3[i - 1] + Y_filler3[i])

X4, Y_filler4 = [], []
with open("BF3vertex.txt") as bf3file:
    ctr = 0
    for line in bf3file:
        if '%' in line:
            ctr += 1
            if ctr != 2:
                continue
            else:
                break
        X4.append(int(line.split(',')[0]))
        Y_filler4.append(int(line.split(',')[2]))    

Y4 = []
for i in range(0, len(Y_filler4)):
    if i == 0:
        Y4.append(Y_filler4[i])
    else:
        Y4.append(Y4[i - 1] + Y_filler4[i])

plt.plot(X, Y, color='blue')
plt.plot(X, X, color='green')
plt.plot(X2, Y2, color='red')
plt.plot(X3, Y3, color = 'purple')
plt.plot(X4, Y4, color='black')
plt.show()
