import matplotlib.pyplot as plt

# dashList = [(5,2),(2,5),(4,10),(3,3,2,2),(5,2,20,2)]

x = [4.53, 16.6, 125.9, 404.2, 466.1, 2413, 3377, 3601, 3601]
y = [0, 1, 2, 3, 4, 5, 6, 7, 10]
plt.step(x, y, label='Acc. Bilevel Decomp.', linestyle='-')

x = [7.9, 16.5, 1230.1, 1911, 2689.5, 3601, 3601]
y = [0, 1, 2, 3, 4, 5, 10]
plt.step(x, y, label='Bilevel Decomp.', linestyle='--', dashes=(5,2))

x = [5.55, 35, 279.9, 2149, 3601, 3601]
y = [0, 1, 2, 3, 4, 10]
plt.step(x, y, label='BARON', linestyle='-.')

x = [7.96, 42.25, 3601, 3601]
y = [0, 1, 2, 10]
plt.step(x, y, label='SCIP', linestyle='--', dashes=(3,3,2,2))

x = [96.19, 147.79, 3601, 3601]
y = [0, 1, 2, 10]
plt.step(x, y, label='ANTIGONE', linestyle=':', color='k')
# plt.xscale('log')

plt.legend()
plt.ylabel("Number of instances solved")
plt.xlabel("Time required to solve instances [s]")
plt.xlim(0, 3600)
plt.ylim(0, 10)

plt.savefig('performance_curve.png')