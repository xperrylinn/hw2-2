import matplotlib.pyplot as plt

# Data for strong scaling
thread_counts_1 = [1, 2, 4, 8, 16, 32, 64]
times_1 = [8.108, 6.941, 6.164, 3.711, 2.514, 2.167, 1.961]

# Data for weak scaling
thread_counts_2 = [1, 2, 4, 8, 16, 32, 64]
particle_counts = [10000, 20000, 40000, 80000, 160000, 320000, 640000]
times_2 = [.514, 1.318, 2.555, 2.954, 4.104, 8.885, 17.869]

fig, axs = plt.subplots(2, 1, figsize=(10, 10))

# Plotting strong scaling
axs[0].plot(thread_counts_1, times_1, marker='o')
axs[0].set_title('Current Commit: Strong Scaling')
axs[0].set_xlabel('Thread Count')
axs[0].set_ylabel('Time (s)')
axs[0].grid(True)

# Plotting weak scaling
axs[1].plot(thread_counts_2, times_2, marker='o')
axs[1].set_title('Current Commit: Weak Scaling')
axs[1].set_xlabel('Thread Count')
axs[1].set_ylabel('Time (s)')

# Annotating the particle count
for i, txt in enumerate(particle_counts):
    axs[1].annotate(f'{txt}', (thread_counts_2[i], times_2[i]))
axs[1].grid(True)

plt.tight_layout()
plt.savefig('scaling_graphs.png')
plt.show()
