import matplotlib.pyplot as plt

x = [4.53, 16.6, 125.9, 404.2, 466.1, 2413, 3377, 3601, 3601]
y = [0, 1, 2, 3, 4, 5, 6, 7, 10]
plt.step(x, y, label='Acc. Bilevel Decomp.', linestyle='-', marker='^', linewidth=3.0, markersize=7)

x = [7.9, 16.5, 1230.1, 1911, 2689.5, 3601, 3601]
y = [0, 1, 2, 3, 4, 5, 10]
plt.step(x, y, label='Bilevel Decomp.', linestyle=':',  marker='s', linewidth=4.0, markersize=6)

x = [5.55, 35, 279.9, 2149, 3601, 3601]
y = [0, 1, 2, 3, 4, 10]
plt.step(x, y, label='BARON', linestyle='-.', marker='x', markersize=6)

x = [7.96, 42.25, 3601, 3601]
y = [0, 1, 2, 10]
plt.step(x, y, label='SCIP', linestyle=':', marker='o', linewidth=2.0, markersize=6)

x = [96.19, 147.79, 3601, 3601]
y = [0, 1, 2, 10]
plt.step(x, y, label='ANTIGONE', linestyle='--', marker='*', markersize=7, color='k')
plt.xscale('log')

plt.legend()
plt.ylabel("Number of instances solved")
plt.xlabel("Time required to solve instances [s]")
plt.xlim(0, 3600)
plt.ylim(0, 10)

plt.savefig('performance_curve.png')

x = [4.53, 16.6, 125.9, 404.2, 466.1, 2413, 3377, 3700, 3700]
y = [0, 1, 2, 3, 4, 5, 6, 7, 10]
plt.step(x, y, label='Acc. Bilevel Decomp.', linestyle='-', marker='^', linewidth=3.0, markersize=7)

x = [7.45, 16.55, 1361.91, 2135.15, 3392.58, 3700, 3700]
y = [0, 1, 2, 3, 4, 5, 10]
plt.step(x, y, label='Acc. Bilevel Decomp. w/o Facility Pruning', linestyle=':',  marker='s', linewidth=4.0, markersize=6)

x = [4.05, 16.94, 31.2, 240.3, 710.28, 2156.03, 3700, 3700]
y = [0, 1, 2, 3, 4, 5, 6, 10]
plt.step(x, y, label='Acc. Bilevel Decomp. w/o Partition Pruning', linestyle='-.', marker='x', markersize=6)

x = [4.21, 16.51, 125.32, 455.26, 530.46, 2439.91, 3700, 3700]
y = [0, 1, 2, 3, 4, 5, 6, 10]
plt.step(x, y, label='Acc. Bilevel Decomp. w/o warm-start', linestyle=':', marker='o', linewidth=2.0, markersize=6)
plt.xscale('log')

plt.legend()
plt.ylabel("Number of instances solved")
plt.xlabel("Time required to solve instances [s]")
plt.xlim(0, 3600)
plt.ylim(0, 10)

plt.savefig('performance_curve_2.png')