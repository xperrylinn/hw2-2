import matplotlib.pyplot as plt

n_values_ten = [1000, 10000, 100000, 200000, 500000, 1000000]
serial_ten = [0.0419525, 0.477981, 7.25282, 20.538, 66.5155, 139.237]
parallel_ten = [0.0538448, 0.273497, 1.90122, 4.27425, 13.8935, 33.8892]

# Plotting
plt.figure(figsize=(10, 6))
plt.loglog(n_values_ten, serial_ten, label='Serial')
plt.loglog(n_values_ten, parallel_ten, label='Parallel')
plt.xlabel('Number of Particles (n)')
plt.ylabel('Runtime (seconds)')
plt.title('Serial + Parallel Simulation runtimes in log-log scale')
plt.legend()
plt.grid(True)
plt.savefig('loglog.png')
plt.show()
