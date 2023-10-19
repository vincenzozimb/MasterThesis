import tensornetwork as tn
import numpy as np


def block(*dimensions):
  '''Construct a new matrix for the MPS with random numbers from 0 to 1'''
  size = tuple([x for x in dimensions])
  return np.random.random_sample(size)

a = tn.Node(block(2,2))
b = tn.Node(block(1,1))
# a = block(2,2)
# b = block(1,2)

c = [a] + [b]

print(a.shape, '\n')
print(b.shape, '\n')

ranks = range(2, 60)