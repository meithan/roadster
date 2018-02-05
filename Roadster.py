# Animation of the Falcon Heavy Roadster's trajectory
import locale
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Orbit import *
from SolarSystem import *
from datetime import datetime, timedelta

plt.style.use('dark_background')

launch_date = datetime(2018,2,6)
dtime = timedelta(seconds=1*86400)

# Roadster's orbit
secs = (launch_date-J2000).total_seconds()
x, y, z, vx, vy, vz = Earth.posvel_at_time(secs)
pos = np.array((x,y,z))
vel = np.array((vx,vy,vz))*1.074
Roadster = CelestialBody("Roadster")
Roadster.orbit = Orbit(primary=Sun)
Roadster.orbit.from_statevectors(pos, vel, secs, J2000)

num_its = 365*4
for it in range(num_its):

  date = launch_date + dtime*it
  secs = (date - J2000).total_seconds()

  fig = plt.figure(figsize=(6,3.375))
  ax = fig.add_subplot(111)
  ax.set_aspect('equal')

  locale.setlocale(locale.LC_ALL, 'en_US')

  Earth.plot_orbit(ax=ax, proj="2D", color="deepskyblue", units=AU)
  Earth.plot_at_time(secs)
  Ex,Ey,Ez,_,_,_ = Earth.posvel_at_time(secs)

  Mars.plot_orbit(ax=ax, proj="2D", color="darkorange", units=AU)
  Mars.plot_at_time(secs)
  Mx,My,Mz,_,_,_ = Mars.posvel_at_time(secs)

  Roadster.plot_orbit(ax=ax, proj="2D", color="red", units=AU)
  Roadster.plot_at_time(secs)
  Rx,Ry,Rz,Rvx,Rvy,Rvz = Roadster.posvel_at_time(secs)

  dist_Earth = sqrt((Rx-Ex)**2 + (Ry-Ey)**2 + (Rz-Ez)**2)
  dist_Mars = sqrt((Rx-Mx)**2 + (Ry-My)**2 + (Rz-Mz)**2)
  dist = sqrt(Rx**2 + Ry**2 + Rz**2)
  speed = sqrt(Rvx**2 + Rvy**2 + Rvz**2)

  plt.xlim(-2,2)
  plt.ylim(-2,2)
  for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(6)
  for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(6)

  plt.text(0.5, 0.95, date.strftime("%d / %b / %Y"), transform=ax.transAxes, ha="center", fontsize=8)

  plt.text(0.03, 0.09, "Distance to ...", transform=ax.transAxes, fontsize=6)
  plt.text(0.03, 0.06, "Earth: %.1f million km" % (dist_Earth/1e9), transform=ax.transAxes, fontsize=6)
  plt.text(0.03, 0.03, "Mars: %.1f million km" % (dist_Mars/1e9), transform=ax.transAxes, fontsize=6)

  plt.text(0.58, 0.06, "Heliocentric distance: %.2f AU" % (dist/AU), transform=ax.transAxes, fontsize=6)
  plt.text(0.58, 0.03, "Heliocentric speed: %.1f km/s" % (speed/1e3), transform=ax.transAxes, fontsize=6)

  plt.subplots_adjust(top=0.99, bottom=0.05, left=0.075, right=0.975, hspace=0.2, wspace=0.2)

  fname = "frames/frame%04i.png" % it
  plt.savefig(fname, dpi=320)
  print("Frame %i/%i" % (it+1, num_its))
  #plt.show()

  plt.close()
