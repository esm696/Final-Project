import numpy as np
import matplotlib.pyplot as plt
phi = (np.sqrt(5.0) - 1.0) / 2.0
epsilon = 1.0
TOLERANCE = 1e-6
MAX_STEPS = 100

V = lambda x: 4.0 * epsilon * (x**12 - x**6)

x = [0.6, None, None, 1.0]
x[1] = x[3] - phi * (x[3] - x[0])
x[2] = x[0] + phi * (x[3] - x[0])

success = False
error = []
for n in range(1, MAX_STEPS + 1):
    f_1 = V(x[1])
    f_2 = V(x[2])
    if f_1 < f_2:
        x[3] = x[2]
        x[2] = x[1]
        x[1] = x[3] - phi * (x[3] - x[0])
    else:
        x[0] = x[1]
        x[1] = x[2]
        x[2] = x[0] + phi * (x[3] - x[0])
        error.append(np.abs((x[3] + x[0]) / 2.0 - 1.0 / 2.0**(1.0 / 6.0)))
        if np.abs(x[3] - x[0]) < TOLERANCE:
            success = True
            break
if success:
    print("Success!")
    print(" t* = %s" % str((x[3] + x[0]) / 2.0))
    print(" f(t*) = %s" % V((x[3] + x[0]) / 2.0))
    print(" number of steps = %s" % n)
else:
    print("Reached maximum number of steps!")

fig = plt.figure()
fig.set_figwidth(fig.get_figwidth() * 2)
axes = fig.add_subplot(1, 2, 1)
axes.loglog(range(1, len(error) + 1), error, 'ro--')
axes.set_title('Convergence to true solution')
axes.set_xlabel("n (step)")
axes.set_ylabel("Absolute Error")
axes = fig.add_subplot(1, 2, 2)
axes.plot((x[3] + x[0]) / 2.0, V((x[3] + x[0]) / 2.0), 'ro')
x = np.linspace(0.0, 1.0, 100)
axes.plot(x, V(x), 'k')
axes.set_title('Potential and Minimum')
axes.set_xlabel("x")
axes.set_ylabel("V(x)")
plt.show()


