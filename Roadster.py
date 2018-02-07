# Animation of the Falcon Heavy Roadster's trajectory
import locale
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Orbit import *
from SolarSystem import *
from datetime import datetime, timedelta
from utils import AU

plt.style.use('dark_background')

launch_date = datetime(2018,2,6)
dtime = timedelta(seconds=1*86400)

# Roadster's orbit
x, y, z, vx, vy, vz = Earth.posvel_at_time(launch_date)
pos = np.array((x,y,z))
vel = np.array((vx,vy,vz))*1.1965
Roadster = CelestialBody("Roadster")
Roadster.orbit = Orbit(primary=Sun)
Roadster.orbit.from_statevectors(pos, vel, launch_date, J2000)

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
print(inner_xs.shape, outer_xs.shape)
xs = np.concatenate((inner_xs, outer_xs))
ys = np.concatenate((inner_ys, outer_ys))

num_its = 365*5
for it in range(num_its):

  date = launch_date + dtime*it

  fig = plt.figure(figsize=(6,3.375))
  ax = fig.add_subplot(111)
  ax.set_aspect('equal')

  locale.setlocale(locale.LC_ALL, 'en_US')

  Earth.plot_orbit(ax=ax, proj="2D", color="deepskyblue", units=AU)
  Earth.plot_at_time(date)
  Ex,Ey,Ez,_,_,_ = Earth.posvel_at_time(date)
  plt.text(Ex/AU, Ey/AU+0.12, "Earth", color="deepskyblue", ha="center", fontsize=6)

  Mars.plot_orbit(ax=ax, proj="2D", color="darkorange", units=AU)
  Mars.plot_at_time(date)
  Mx,My,Mz,_,_,_ = Mars.posvel_at_time(date)
  plt.text(Mx/AU, My/AU-0.20,"Mars", color="darkorange", ha="center", fontsize=6)

  Roadster.plot_orbit(ax=ax, proj="2D", color="red", units=AU)
  Roadster.plot_at_time(date)
  Rx,Ry,Rz,Rvx,Rvy,Rvz = Roadster.posvel_at_time(date)
  plt.text(Rx/AU, Ry/AU-0.20, "Roadster", color="red", ha="center", fontsize=6)

  dist_Earth = sqrt((Rx-Ex)**2 + (Ry-Ey)**2 + (Rz-Ez)**2)
  dist_Mars = sqrt((Rx-Mx)**2 + (Ry-My)**2 + (Rz-Mz)**2)
  dist = sqrt(Rx**2 + Ry**2 + Rz**2)
  speed = sqrt(Rvx**2 + Rvy**2 + Rvz**2)

  # Plot main asteroid belt
  ax.fill(xs, ys, color='0.075', lw=0)
  plt.text(-1.9, -1.9, "Asteroid belt", color='0.2', fontsize=6, ha="center")

  ax.plot([0], [0], "x", ms=3, color="0.2")
  plt.text(0, 0.1, "Sun", color="0.2", fontsize=6, ha="center")

  plt.xlim(-3,3)
  plt.ylim(-3.5,2.5)
  for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(6)
  for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(6)

  plt.text(0.5, 0.95, date.strftime("%d / %b / %Y"), transform=ax.transAxes, ha="center", fontsize=8)

  plt.text(0.03, 0.09, "Distance to ...", transform=ax.transAxes, fontsize=6)
  plt.text(0.03, 0.06, "Earth: %.1f million km" % (dist_Earth/1e9), transform=ax.transAxes, fontsize=6)
  plt.text(0.03, 0.03, "Mars: %.1f million km" % (dist_Mars/1e9), transform=ax.transAxes, fontsize=6)

  #plt.text(0.58, 0.09, "Roadster's ...", transform=ax.transAxes, fontsize=6)
  plt.text(0.58, 0.06, "Heliocentric distance: %.2f AU" % (dist/AU), transform=ax.transAxes, fontsize=6)
  plt.text(0.58, 0.03, "Heliocentric speed: %.1f km/s" % (speed/1e3), transform=ax.transAxes, fontsize=6)

  plt.subplots_adjust(top=0.99, bottom=0.05, left=0.075, right=0.975, hspace=0.2, wspace=0.2)

  fname = "frames/frame%04i.png" % it
  plt.savefig(fname, dpi=320)
  print("Frame %i/%i" % (it+1, num_its))
  #plt.show()

  plt.close()
