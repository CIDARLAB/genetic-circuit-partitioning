import matplotlib.pyplot as plt 
import math

Tmax = 21  # starting temperature
C = 1e-4
TiList, PList = [], []
steps = 10000
for i in range(steps):
	Ti = Tmax*(math.exp(-C*i))
	P = math.exp(-1/Ti)
	TiList.append(Ti)
	PList.append(P)

Plist_test = []
xlist = []
i = 3
Ti = Tmax*(math.exp(-C*i))
P = math.exp(-1/Ti)
print(Ti)
print(P)

# for deltaT in range(-5, 6):
# 	xlist.append(deltaT)
# 	Ti = Tmax*(math.exp(-C*i))
# 	print('Ti', Ti)
# 	Plist_test.append(math.exp(deltaT/Ti))
# # plt.plot(range(steps), [math.exp(i) for i in range(steps)])
# # plt.plot(range(steps), TiList)
# # plt.plot(range(steps), PList)

# for t in range(1000):
# 	Ti = Tmax*(math.exp(-C*t))
# 	P = math.exp(-1/Ti)
# 	Plist_test.append(P)


plt.plot(range(1000), Plist_test)
plt.show()
