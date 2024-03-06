import matplotlib.pyplot as plt

# Thread counts
threads = [1, 2, 4, 8, 16, 32, 64]

# Serial runtime
serial_runtime = 7.3082

# Parallel runtimes for different thread counts
parallel_runtimes = {
    1: 8.09461,
    2: 6.92824,
    4: 6.29216,
    8: 3.7479,
    16: 2.50816,
    32: 2.19544,
    64: 1.99988
}

# Calculate speedup
ideal_speedup = [parallel_runtimes[1] / p for p in threads]  # Ideal speedup where speedup = T1 / p

speedup = [parallel_runtimes[p] / p for p in threads]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(threads, speedup, marker='o', linestyle='-', color='b', label='Speedup')
plt.plot(threads, ideal_speedup, linestyle='--', color='r', label='Ideal Speedup')
plt.xlabel('Number of Threads')
plt.ylabel('Speedup')
plt.title('Speedup Plot')
plt.legend()
plt.grid(True)
plt.xticks(threads)
plt.savefig('speedup.png')
plt.show()
