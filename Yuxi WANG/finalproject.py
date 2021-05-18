import numpy as np
import matplotlib.pyplot as plt

# use golden section method
phi = (np.sqrt(5) - 1) / 2  # phi=golden ratio 0.6180
epsilon = 1.57  # depth of potential well
TOLERANCE = 1e-6
MAX_STEPS = 100

V = lambda x: 4 * epsilon * (x ** 12 - x ** 6)

x = [0.6, None, None, 1.0]
x[1] = x[3] - phi * (x[3] - x[0])
x[2] = x[0] + phi * (x[3] - x[0])

success = False
error = []
for n in range(1, MAX_STEPS + 1):
    f_1 = V(x[1])
    f_2 = V(x[2])

    if f_1 < f_2:  # the minimum is between  x[0]  and  x[2]
        x[3] = x[2]
        x[2] = x[1]
        x[1] = x[3] - phi * (x[3] - x[0])

    else:  # the minimum is between  x[1]  and  x[3]
        x[0] = x[1]
        x[1] = x[2]
        x[2] = x[0] + phi * (x[3] - x[0])

        error.append(np.abs((x[3] + x[0]) / 2 - 1 / 2 ** (1 / 6)))
        if np.abs(x[3] - x[0]) < TOLERANCE:
            success = True
            break

if success:
    print("Success")
    print(" t* = %s" % str((x[3] + x[0]) / 2))  # t*=x with minimum potential
    print(" f(t*) = %s" % V((x[3] + x[0]) / 2))  # f(t*)=minimum potential
    print(" number of steps = %s" % n)
else:
    print("Reached maximum number of steps")

fig = plt.figure()
fig.set_figwidth(fig.get_figwidth() * 2)
axes = fig.add_subplot(1, 2, 1)
axes.loglog(range(1, len(error) + 1), error, 'ro--')
axes.set_title('Convergence to true solution')
axes.set_xlabel("n (step)")
axes.set_ylabel("Absolute Error")
axes = fig.add_subplot(1, 2, 2)
axes.plot((x[3] + x[0]) / 2, V((x[3] + x[0]) / 2), 'ro')  # illustrate the minima point
x = np.linspace(0, 1.0, 100)
axes.plot(x, V(x), 'k')
axes.set_title('Potential and Minimum')
axes.set_xlabel("x")
axes.set_ylabel("V(x)")
plt.show()


