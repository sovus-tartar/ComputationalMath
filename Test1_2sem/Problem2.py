import numpy as np
import matplotlib.pyplot as plt

# Define the complex plane
x = np.linspace(-5, 5, 500)
y = np.linspace(-5, 5, 500)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y

# Define the complex set
C1 = (np.abs(1 + 0.81 * Z) <= np.abs(1 - 0.19 * Z)) #or (np.abs(X) <= 0.05) or (np.abs(Y) <= 0.05)
# C2 = (np.imag(Z) > 0)
# Plot the complex set
plt.figure(figsize=(6,6))
plt.contourf(X, Y, C1, cmap='coolwarm')
# plt.contourf(X, Y, C2, cmap='coolwarm')
plt.colorbar()
plt.show()