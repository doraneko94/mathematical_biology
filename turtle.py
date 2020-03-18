import numpy as np
from matplotlib import pyplot as plt

L = np.array([[   0.0,    0.0,    0.0,    0.0,  127.0,    4.0,   80.0],
              [0.6747, 0.7370,    0.0,    0.0,    0.0,    0.0,    0.0],
              [   0.0, 0.0486, 0.6610,    0.0,    0.0,    0.0,    0.0],
              [   0.0,    0.0, 0.0147, 0.6907,    0.0,    0.0,    0.0],
              [   0.0,    0.0,    0.0, 0.0518,    0.0,    0.0,    0.0],
              [   0.0,    0.0,    0.0,    0.0, 0.8091,    0.0,    0.0],
              [   0.0,    0.0,    0.0,    0.0,    0.0, 0.8091, 0.8089]])

l, V = np.linalg.eig(L)
l1 = np.abs(l).max()
V1 = np.real(V[:, np.abs(l).argmax()])
v1 = V1.sum()
stable = V1 / v1

print("==========================================================")
print("Dominant Eigenvalue in a Narrow Sense: lambda_1 = {:.5f}".format(l1))
print("")
print("Stable Stage Distribution: V1/v1")
print("(eggs, hatchlings)  = {:.2f}%".format(stable[0] * 100))
print("(small juveniles)   = {:.2f}%".format(stable[1] * 100))
print("(large juveniles)   = {:.2f}%".format(stable[2] * 100))
print("(subadults)         = {:.2f}%".format(stable[3] * 100))
print("(novice breeders)   = {:.2f}%".format(stable[4] * 100))
print("(1st-yr remigrants) = {:.2f}%".format(stable[5] * 100))
print("(mature breeders)   = {:.2f}%".format(stable[6] * 100))
print("==========================================================")

num = 1000000
x = np.arange(51)

turtle_v = stable.reshape(-1, 1) * num
turtle_v = turtle_v.astype(np.int32)
y1 = [turtle_v.sum()]
for _ in range(50):
    turtle_v = np.dot(L, turtle_v).astype(np.int32)
    y1.append(turtle_v.sum())
y1 = np.array(y1)

turtle_v = stable.reshape(-1, 1) * num
turtle_v = turtle_v.astype(np.int32)
y2 = [turtle_v.sum()]
for _ in range(50):
    turtle_v = (turtle_v * l1).astype(np.int32)
    y2.append(turtle_v.sum())
y2 = np.array(y2)

plt.plot(x, y1)
plt.plot(x, y2)
plt.legend(["simulation", "approximation"])
plt.xlabel("Year")
plt.ylabel("Total Population")
plt.title("Crouse et al. (1987)")
plt.tight_layout()
plt.savefig("turtle.png")