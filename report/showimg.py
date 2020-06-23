import numpy as np
import matplotlib.pyplot as plt


serial1 = np.array([0.00193037, 0.005019587, 0.015565,
                    0.0255725, 0.110831, 0.157498])
parallel1 = np.array([0.009078, 0.0239482, 0.0709564,
                      0.364968, 0.357466, 0.734400])
parallel32 = np.array([0.0119314, 0.0136889, 0.00907152,
                       0.0165627, 0.0437735, 0.0436356])

bar_width = 0.25
x = np.arange(6)
xlabel = ['case7', 'case8', 'case9', 'case10', 'case11', 'case12']

speedup = serial1/parallel32


plt.bar(x, serial1, bar_width, label="serial1")
#plt.bar(x, parallel1, bar_width, label="parallel1")
plt.bar(x + bar_width, parallel32, bar_width, label="parallel32")
plt.xticks(x + bar_width / 2, xlabel)
plt.xlabel('Case')
plt.ylabel('Time(s)')
plt.title('parallel32-serial1 Compare')
plt.legend()
ax = plt.twinx()
ax.plot(x + bar_width / 2, speedup, c="black")
ax.set_ylabel("SpeedUp")
plt.savefig("./res1.png")
# plt.show()
