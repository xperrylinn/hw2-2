import matplotlib.pyplot as plt

# Data for strong scaling
thread_counts = [1, 2, 4, 8, 16, 32, 64]
times_1 = [3.26996, 2.37017, 1.4263, 0.837236, 0.460013, 0.332347, 0.333075]

# Data for weak scaling
#thread_counts_2 = [1, 2, 4, 8, 16, 32, 64]
particle_counts = [10000, 20000, 40000, 80000, 160000, 320000, 640000]
times_2 = [3.26996, 4.9118, 6.397, 8.18432, 16.2359, 19.223, 29.1091]

fig, axs = plt.subplots(2, 1, figsize=(10, 10))

# Plotting strong scaling
axs[0].plot(thread_counts, times_1, marker='o')
axs[0].set_title('Strong Scaling (1 Node)')
axs[0].set_xlabel('Thread Count')
axs[0].set_ylabel('Time (s)')
axs[0].grid(True)

# Plotting weak scaling
axs[1].plot(thread_counts, times_2, marker='o')
axs[1].set_title('Weak Scaling (1 Node)')
axs[1].set_xlabel('Thread Count')
axs[1].set_ylabel('Time (s)')

# Annotating the particle count
for i, txt in enumerate(particle_counts):
    axs[1].annotate(f'{txt}', (thread_counts[i], times_2[i]))
axs[1].grid(True)

plt.tight_layout()
plt.savefig('scaling_graphs.png')
plt.show()
