import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# Set up figure
plt.figure(figsize=(6, 6)) 
plt.gca().set_aspect('equal', adjustable='datalim') 

# Define custom color list
colors = ['red', 'blue', 'green', 'orange', 'magenta']




# Load and plot data with specified colors
log = np.genfromtxt("Data/Error/re1000/error_re=1000_linear.csv", delimiter=',')
dH = log[:, 0]
error = log[:, 1]
plt.plot(np.log2(dH), np.log2(error), marker='.', color=colors[0], label='Δt=Δh - Rₑ = 1000')

log = np.genfromtxt("Data/Error/re1000/error_re=1000_quadratic.csv", delimiter=',')
dH = log[:, 0]
error = log[:, 1]
plt.plot(np.log2(dH), np.log2(error), marker='.', color=colors[1], label='Δt=Δh² - Rₑ = 1000')

log = np.genfromtxt("Data/Error/re10000/error_re=10000_linear.csv", delimiter=',')
dH = log[:, 0]
error = log[:, 1]
plt.plot(np.log2(dH), np.log2(error), marker='.', color=colors[2], label='Δt=Δh - Rₑ = 10000')

log = np.genfromtxt("Data/Error/re10000/error_re=10000_quadratic.csv", delimiter=',')
dH = log[:, 0]
error = log[:, 1]
plt.plot(np.log2(dH), np.log2(error), marker='.', color=colors[3], label='Δt=Δh² - Rₑ = 10000')

plt.grid(True)
#plt.title("Logarithmic Error Decay in Taylor-Green Vortex", fontsize=16, fontweight='bold')
plt.xlabel("Log2(N)")
plt.ylabel("Log2(Max Absolute Error)")
plt.xlim(2, 6)
plt.ylim(-8, -2)

# Annotate convergence slope
x_annotation = 4.0
y_annotation = -12
plt.annotate("k = 2", xy=(x_annotation-0.5, y_annotation+0.5),
             xytext=(x_annotation-0.4, y_annotation-0.35))

plt.annotate("", xy=(x_annotation, y_annotation),
             xytext=(x_annotation-0.5, y_annotation+1),
             arrowprops=dict(facecolor='black', arrowstyle='->'))
plt.annotate("", xy=(x_annotation, y_annotation),
             xytext=(x_annotation-0.5, y_annotation),
             arrowprops=dict(facecolor='black', arrowstyle='-'))
plt.annotate("", xy=(x_annotation-0.5, y_annotation),
             xytext=(x_annotation-0.5, y_annotation+1),
             arrowprops=dict(facecolor='black', arrowstyle='-'))

plt.plot(x_annotation - 0.5, y_annotation + 1, '.', color='black')
plt.plot(x_annotation - 0.5, y_annotation, '.', color='black')
plt.plot(x_annotation, y_annotation, '.', color='black')

plt.legend(loc='best', framealpha=0.3)
plt.savefig("error_log.png")
plt.show()
plt.close()
