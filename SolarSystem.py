# Celestial bodies of the solar system
from CelestialBody import CelestialBody
from Orbit import Orbit
from utils import AU, J2000

# mu given in MKS
Sun = CelestialBody("Sun", mu=1.32712440018e20)
Sun.radius = 695700e3

# Solar system orbital elements
# The following table was obtained from:
# http://www.met.rdg.ac.uk/~ross/Astronomy/Planets.html
# It lists: a (AU), e, i, Omega/LAN, ~omega, L0
# where ~omega is the longitude of the periapsis and L0 is the mean longitude
# at epoch. The argument of periapsis omega and mean anomaly at epoch M0 can
# then be obtained as:
#   omega = ~omega - LAN
#   M0 = L0 - ~omega
# The epoch is J2000.
elements = [
["Mercury", 0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084],
["Venus", 0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973],
["Earth", 1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435],
["Mars", 1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332],
["Jupiter", 5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438],
["Saturn", 9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432],
["Uranus", 19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218],
["Neptune", 30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003],
["Pluto", 39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881]
]

ss_planets = {}
for name,a,e,i,LAN,long_peri,L0 in elements:

  planet = CelestialBody(name)
  planet.orbit = Orbit(primary=Sun)
  arg_peri = long_peri - LAN
  M0 = L0 - long_peri
  planet.orbit.from_elements(a=a*AU, e=e, i=i, arg=arg_peri, LAN=LAN, M=M0, time=J2000)

  ss_planets[name] = planet

ss_planets["Earth"].mu = 3.986004418e14
ss_planets["Earth"].radius = 6370.0e3

Mercury = ss_planets["Mercury"]
Venus = ss_planets["Venus"]
Earth = ss_planets["Earth"]
Mars = ss_planets["Mars"]
Jupiter = ss_planets["Jupiter"]
Saturn = ss_planets["Saturn"]
Uranus = ss_planets["Uranus"]
Neptune = ss_planets["Neptune"]
Pluto = ss_planets["Pluto"]

# ==============================================

# Plots the current positions of the planets
if __name__ == "__main__":

  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import Axes3D
  plt.style.use('dark_background')
  from datetime import datetime

  time = datetime.utcnow()
  # time = datetime(2018,1,1)

  proj = "3D"   # "3D" or "2D"
  inner_only = True   # If True, plots only the inner planets

  fig = plt.figure(figsize=(6,6))
  if proj == "3D":
    ax = fig.add_subplot(111, projection='3d')
  elif proj == "2D":
    ax = fig.add_subplot(111)
  ax.set_aspect('equal')

  Mercury.plot_orbit(ax=ax, proj=proj, color="darkgray", show_axes=True, units=AU)
  Mercury.plot_at_time(time)
  Venus.plot_orbit(ax=ax, proj=proj, color="peru", units=AU)
  Venus.plot_at_time(time)
  Earth.plot_orbit(ax=ax, proj=proj, color="darkcyan", units=AU)
  Earth.plot_at_time(time)
  Mars.plot_orbit(ax=ax, proj=proj, color="firebrick", units=AU)
  Mars.plot_at_time(time)

  if not inner_only:
    Jupiter.plot_orbit(ax=ax, proj=proj, color="chocolate", units=AU)
    Jupiter.plot_at_time(time)
    Saturn.plot_orbit(ax=ax, proj=proj, color="darksalmon", units=AU)
    Saturn.plot_at_time(time)
    Uranus.plot_orbit(ax=ax, proj=proj, color="steelblue", units=AU)
    Uranus.plot_at_time(time)
    Neptune.plot_orbit(ax=ax, proj=proj, color="royalblue", units=AU)
    Neptune.plot_at_time(time)
    Pluto.plot_orbit(ax=ax, proj=proj, color="gray", units=AU)
    Pluto.plot_at_time(time)
    plt.xlim(-48,48)
    plt.ylim(-48,48)

  # for name in ss_planets:
  #   planet = ss_planets[name]
  #   print("%s %.3f x %.3f AU, %.1f°" % (name.ljust(8), planet.orbit.rpe/AU, planet.orbit.rap/AU, planet.orbit.i))

  if proj == "2D":
    plt.text(0.97, 0.95, now.strftime("%Y-%m-%d"), transform=ax.transAxes, ha="right")

  ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
  ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
  ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
  ax.grid(False)

  plt.subplots_adjust(top=0.99, bottom=0.05, left=0.075, right=0.975, hspace=0.2, wspace=0.2)

  plt.show()
