# import necessary packages
import numpy as np
import tensornetwork as tn

# create a random tensor
def block(*dimensions):
  '''Construct a tensor with random numbers from 0 to 1'''
  size = tuple([x for x in dimensions])
  return np.random.random_sample(size)

