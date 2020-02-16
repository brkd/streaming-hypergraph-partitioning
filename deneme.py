import matplotlib.pyplot as plt

X, Y_filler = [], []
with open("N2Pvertex.txt") as n2pfile:
    for line in n2pfile:
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
    for line in bffile:
        X2.append(int(line.split(',')[0]))
        Y_filler2.append(int(line.split(',')[2]))    

Y2 = []
for i in range(0, len(Y_filler2)):
    if i == 0:
        Y2.append(Y_filler2[i])
    else:
        Y2.append(Y2[i - 1] + Y_filler2[i])

plt.plot(X, Y, color='blue')
plt.plot(X, X, color='green')
plt.plot(X2, Y2, color='red')
plt.show()

    
