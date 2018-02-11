# Animation of the Falcon Heavy Roadster's trajectory
# Reads heliocentric positions from datafiles obtained from JPL Horizons
from datetime import datetime, timedelta
import locale
from math import pi, sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# ==============================================================================
# CONFIGURATION

# Data files with JPL Horizons heliocentric positions
Mercury_fname = "horizons_Mercury.txt"
Venus_fname = "horizons_Venus.txt"
Earth_fname = "horizons_Earth.txt"
Mars_fname = "horizons_Mars.txt"
Roadster_fname = "horizons_Roadster.txt"

# Some special dates
launch_date = datetime(2018, 2, 6)
first_crossing = datetime(2018, 7, 12)
first_aphelion = datetime(2018, 11, 9)
second_crossing = datetime(2019, 3, 28)
first_perihelion = datetime(2019, 8, 15)
closest_Mars_approach = datetime(2020, 10, 7)
closest_Earth_approach = datetime(2021, 3, 29)

# Orbital periods
Mercury_P = 176
Venus_P = 243
Earth_P = 365
Mars_P = 687
Roadster_P = 557

# Directory to output movie frames
frames_dir = "frames/"

# Colors
Mercury_color = "gray"
Venus_color = "tan"
Earth_color = "dodgerblue"
Mars_color = "darkorange"
Roadster_color = "#ff0021"   # ff0021 = Bright red
Sun_color = "gold"

# Default plot params
fontsize = 18
msize = 60

# ==============================================================================

AU_km = 149597870.7

# Load positions from data files
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
      X = float(line[4:26])
      Y = float(line[30:52])
      Z = float(line[56:78])
      line = f.readline()
      VX = float(line[4:26])
      VY = float(line[30:52])
      VZ = float(line[56:78])
      data.append((date,(X/AU_km,Y/AU_km,Z/AU_km),(VX,VY,VZ)))
  return zip(*data)

dates, Roadster_XYZ, Roadster_V = read_JPL_data(Roadster_fname)
Roadster_XYZ = np.array(Roadster_XYZ)
Mercury_XYZ = np.array(list(read_JPL_data(Mercury_fname))[1])
Venus_XYZ = np.array(list(read_JPL_data(Venus_fname))[1])
Earth_XYZ = np.array(list(read_JPL_data(Earth_fname))[1])
Mars_XYZ = np.array(list(read_JPL_data(Mars_fname))[1])

# Asteroid belt
belt_min = 2.06
belt_max = 3.27
numpts = 200
theta = np.linspace(0, 2*pi, numpts, endpoint=True)
inner_xs = belt_min*np.cos(theta)
inner_ys = belt_min*np.sin(theta)
outer_xs = belt_max*np.cos(theta)
outer_ys = belt_max*np.sin(theta)
outer_xs[:] = outer_xs[::-1]
outer_ys[:] = outer_ys[::-1]
xs = np.concatenate((inner_xs, outer_xs))
ys = np.concatenate((inner_ys, outer_ys))

plt.style.use('dark_background')

if not frames_dir.endswith("/"):
  frames_dir += "/"
bbox = dict(pad=0, fc=(0,0,0,0.7), ec="none", lw=0)

num_frames = len(dates) + 5
frame = 0
i = 0
#i = 1140
while i < len(dates):

  date = dates[i]
  Mex, Mey, Mez = Mercury_XYZ[i]
  Vx, Vy, Vz = Venus_XYZ[i]
  Ex, Ey, Ez = Earth_XYZ[i]
  Mx, My, Mz = Mars_XYZ[i]
  Rx, Ry, Rz = Roadster_XYZ[i]
  Rvx, Rvy, Rvz = Roadster_V[i]

  fig = plt.figure(figsize=(11.3,10.8))
  ax = fig.add_subplot(111)
  ax.set_aspect('equal')
  locale.setlocale(locale.LC_ALL, 'en_US')

  # Mercury
  plt.scatter([Mex], [Mey], s=msize, marker="o", color=Mercury_color, zorder=10)
  i1 = i
  i2 = i + Mercury_P + 1
  if i2 >= len(dates):
    i2 = len(dates) - 1
    i1 = i2 - Mercury_P - 1
  plt.plot(Mercury_XYZ[i1:i2,0], Mercury_XYZ[i1:i2,1], color=Mercury_color)
  plt.text(Mex, Mey+0.08, "Mercury", color=Mercury_color, ha="center", fontsize=fontsize, bbox=bbox)

  # Venus
  plt.scatter([Vx], [Vy], s=msize, marker="o", color=Venus_color, zorder=10)
  i1 = i
  i2 = i + Venus_P + 1
  if i2 >= len(dates):
    i2 = len(dates) - 1
    i1 = i2 - Venus_P - 1
  plt.plot(Venus_XYZ[i1:i2,0], Venus_XYZ[i1:i2,1], color=Venus_color)
  plt.text(Vx, Vy+0.08, "Venus", color=Venus_color, ha="center", fontsize=fontsize, bbox=bbox)

  # Earth
  plt.scatter([Ex], [Ey], s=msize, marker="o", color=Earth_color, zorder=10)
  i1 = i
  i2 = i + Earth_P + 1
  if i2 >= len(dates):
    i2 = len(dates) - 1
    i1 = i2 - Earth_P - 1
  plt.plot(Earth_XYZ[i1:i2,0], Earth_XYZ[i1:i2,1], color=Earth_color)
  plt.text(Ex, Ey+0.08, "Earth", color=Earth_color, ha="center", fontsize=fontsize+2, bbox=bbox)

  # Mars
  plt.scatter([Mx], [My], s=msize, marker="o", color=Mars_color, zorder=10)
  i1 = i
  i2 = i + Mars_P + 1
  if i2 >= len(dates):
    i2 = len(dates) - 1
    i1 = i2 - Mars_P - 1
  plt.plot(Mars_XYZ[i1:i2,0], Mars_XYZ[i1:i2,1], color=Mars_color)
  plt.text(Mx, My+0.08, "Mars", color=Mars_color, ha="center", fontsize=fontsize+2, bbox=bbox)

  # Roadster
  plt.scatter([Rx], [Ry], s=msize, marker="o", color=Roadster_color, zorder=10)
  i1 = i
  i2 = i + Roadster_P + 1
  if i2 >= len(dates):
    i2 = len(dates) - 1
    i1 = i2 - Roadster_P - 1
  plt.plot(Roadster_XYZ[i1:i2,0], Roadster_XYZ[i1:i2,1], color=Roadster_color, lw=3)
  plt.text(Rx, Ry-0.12, "Roadster", color=Roadster_color, ha="center", fontsize=fontsize+2, bbox=bbox)

  # Sun
  ax.scatter([0], [0], marker="o", s=msize, color="gold")
  plt.text(0, 0.08, "Sun", color="gold", fontsize=fontsize-2, ha="center")


  dist_Earth = sqrt((Rx-Ex)**2+(Ry-Ey)**2+(Rz-Ez)**2)
  dist_Mars = sqrt((Rx-Mx)**2+(Ry-My)**2+(Rz-Mz)**2)
  helio_dist = sqrt(Rx**2+Ry**2+Rz**2)
  helio_speed = sqrt(Rvx**2+Rvy**2+Rvz**2)

  # Plot main asteroid belt
  # ax.fill(xs, ys, color='0.05', lw=0)
  # plt.text(-1.9, -1.9, "Asteroid belt", color='0.2', fontsize=6, ha="center")

  if frame != 0:
    plt.text(0.5, 0.95, date.strftime("%d / %b / %Y"), transform=ax.transAxes, ha="center", fontsize=fontsize+4)

  plt.text(0.02, 0.08, "Distance to ...", transform=ax.transAxes, fontsize=fontsize)
  plt.text(0.02, 0.05, "Earth: %.1f million km" % (dist_Earth*AU_km/1e6), transform=ax.transAxes, fontsize=fontsize)
  plt.text(0.02, 0.02, "Mars: %.1f million km" % (dist_Mars*AU_km/1e6), transform=ax.transAxes, fontsize=fontsize)

  #plt.text(0.58, 0.09, "Roadster's ...", transform=ax.transAxes, fontsize=6)
  plt.text(0.58, 0.05, "Heliocentric distance: %.3f AU" % (helio_dist), transform=ax.transAxes, fontsize=fontsize)
  plt.text(0.58, 0.02, "Heliocentric speed: %.2f km/s" % (helio_speed), transform=ax.transAxes, fontsize=fontsize)

  if frame == 0:
    plt.annotate("Launch\n%s" % launch_date.strftime("%b %d, %Y"), xy=(Ex, Ey), xytext=(-130,40), textcoords="offset pixels", fontsize=fontsize, bbox=bbox, ha="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
  elif date == first_crossing:
    plt.annotate("Mars orbit crossing\n%s" % date.strftime("%b %d, %Y"), xy=(Rx, Ry), xytext=(-80,80), textcoords="offset pixels", fontsize=fontsize, bbox=bbox, ha="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
  elif date == first_aphelion:
    plt.annotate("Aphelion, 1.66 AU\n%s" % date.strftime("%b %d, %Y"), xy=(Rx, Ry), xytext=(-80,80), textcoords="offset pixels", fontsize=fontsize, bbox=bbox, ha="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
  elif date == second_crossing:
    plt.annotate("Mars orbit crossing\n%s" % date.strftime("%b %d, %Y"), xy=(Rx, Ry), xytext=(10,80), textcoords="offset pixels", fontsize=fontsize, bbox=bbox, ha="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
  elif date == first_perihelion:
    plt.annotate("Perihelion, 0.99 AU\n%s" % date.strftime("%b %d, %Y"), xy=(Rx, Ry), xytext=(10,80), textcoords="offset pixels", fontsize=fontsize, bbox=bbox, ha="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
  elif date == closest_Mars_approach:
    plt.annotate("Closest Mars approach\n7.4 million km\n%s" % date.strftime("%b %d, %Y"), xy=(Mx, My), xytext=(10,80), textcoords="offset pixels", fontsize=fontsize, bbox=bbox, ha="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
  elif date == closest_Earth_approach:
    plt.annotate("Closest Earth approach\n41.1 million km\n%s" % date.strftime("%b %d, %Y"), xy=(Ex, Ey), xytext=(40,120), textcoords="offset pixels", fontsize=fontsize, bbox=bbox, ha="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))


  plt.xlim(-2,2)
  plt.ylim(-2,2)
  for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(fontsize-2)
  for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(fontsize-2)
  plt.xlabel("AU", fontsize=fontsize-2)
  plt.ylabel("AU", fontsize=fontsize-2)

  plt.subplots_adjust(top=0.98, bottom=0.06, left=0.07, right=0.99, hspace=0.2, wspace=0.19)

  fname = frames_dir + "frame%04i.png" % frame
  #plt.savefig(fname, dpi=150)
  plt.savefig(fname)

  print("Frame %i/%i" % (frame, num_frames-1))
  #plt.show()

  plt.close()

  if frame != 0:
    i += 1
  frame += 1
