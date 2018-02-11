# Plots the Roadster's distances to the Sun, the Earth and Mars vs time
import locale
from math import sqrt
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

locale.setlocale(locale.LC_ALL, 'en_US.utf8')

def read_JPL_data(fname):
  start = False
  f = open(fname)
  data = []
  while True:
    line = f.readline()
    if "$$SOE" in line:
      start = True
    elif "$$EOE" in line:
      break
    elif start:
      date = datetime.strptime(line[25:36], "%Y-%b-%d")
      line = f.readline()
      if not line: break
      X = float(line[4:26])
      Y = float(line[30:52])
      Z = float(line[56:78])
      data.append((date,(X/1e6,Y/1e6,Z/1e6),sqrt(X**2+Y**2+Z**2)/1e6))
      line = f.readline()
  return zip(*data)

dates, Roadster_pos, Roadster_helio = read_JPL_data("horizons_Roadster.txt")
_, Earth_pos,_ = read_JPL_data("horizons_Earth.txt")
_, Mars_pos,_ = read_JPL_data("horizons_Mars.txt")
Mars_dist = []
Earth_dist = []
Mars_xy = []
Mars_z = []
for i in range(len(dates)):
  x,y,z = Roadster_pos[i]
  Ex,Ey,Ez = Earth_pos[i]
  Mx,My,Mz = Mars_pos[i]
  Mars_dist.append(sqrt((x-Mx)**2+(y-My)**2+(z-Mz)**2))
  Mars_z.append(abs(z-Mz))
  Mars_xy.append(sqrt((x-Mx)**2+(y-My)**2))
  Earth_dist.append(sqrt((x-Ex)**2+(y-Ey)**2+(z-Ez)**2))

plt.figure(figsize=(12,4))
plt.title("Roadster's distance to ...")
plt.plot(dates, Roadster_helio, "-", color="k", label="Sun")
plt.plot(dates, Earth_dist, "-", color="b", label="Earth")
plt.plot(dates, Mars_dist, "-", color="r", label="Mars")
plt.plot(dates, Mars_xy, "--", color="r", label="Mars xy")
plt.plot(dates, Mars_z, ":", color="r", label="Mars z")

plt.axhline(0,ls=":",color="k")
plt.legend(loc="upper left")
plt.ylabel("Million kilometers")

plt.subplots_adjust(top=0.9,
bottom=0.095,
left=0.055,
right=0.98,
hspace=0.2,
wspace=0.2)

plt.show()
