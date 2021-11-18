import timeit
from copy import deepcopy
# import genetic_partition_test as gp 
import networkx as nx
import numpy as np
import ujson

m1 = np.zeros(shape=(6,6))
m2 = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])

for i in range(m2.shape[0]):
	for j in range(m2.shape[1]):
		m1[i][j] = m2[i][j]
print(m1)