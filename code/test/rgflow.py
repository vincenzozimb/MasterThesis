# import libraries
import numpy as np
import matplotlib.pyplot as plt
import os

print("Hello")

# functions
def f(x, y):
    return 0.25 * np.log(np.cosh(2*x+y) * np.cosh(2*x-y) / (np.cosh(y)**2))
                  
def g(x, y):
    return y * 0.5 * np.log(np.cosh(2*x+y) / np.cosh(2*x-y))

def inverse(x):
    return 0.5 * np.log((1 + x) / (1 - x))

# Create a grid of coordinates
x = np.linspace(0, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)

# Define vector field components (here, a rotational field)
u = f(X, Y)
v = g(X, Y)


# Create a streamplot
plt.streamplot(X, Y, np.tanh(u), np.tanh(v), density=1.5, linewidth=1, color='blue')
# plt.streamplot(np.tanh(X), np.tanh(Y), u, v, density=1.5, linewidth=1, color='blue')

# Customize the plot as needed
plt.title('Streamplot Example')
plt.xlabel('$K_1$')
plt.ylabel('$K_2$')
plt.grid(True)
plt.show()