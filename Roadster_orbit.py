# Annotated diagram of the Roadster's orbit
import locale
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Orbit import *
from SolarSystem import *

# ==============================================================================
# CONFIGURATION

# Figure filename
fname = "Roadster_orbit.png"

# Dates
# Crossing dates are approximate
launch_date = datetime(2018, 2, 6)
crossing1_date = datetime(2018, 7, 12)
crossing2_date = datetime(2019, 3, 28)
diagram_date = launch_date

# Aphelion and perihelion distances and dates
# These are obtained directly from JPL Horizons heliocentric positions
# instead of osculating elements (as those change over time)
aphelion_date = datetime(2018, 11, 9)
aphelion_dist = 1.66   # AU
perihelion_date = datetime(2019, 8, 15)
perihelion_dist = 0.99   # AU

# Colors
Mercury_color = "gray"
Venus_color = "tan"
Earth_color = "dodgerblue"
Mars_color = "darkorange"
Roadster_color = "#ff0021"   # ff0021 = Bright red
Sun_color = "gold"

# Roadster's orbit
# JPL Horizons heliocentric osculating elements for 2018-May-01 (obtained 2018-Feb-10)
# Elements are taken in the future to account for residual Earth perturbation
# Only A, EC, IN, W, OM, MA and the time of the elset are used
EC = 2.557131166888885E-01
QR = 9.861410125199538E-01
IN = 1.077436770097831E+00
OM = 3.171104011083243E+02
W  = 1.774800604644223E+02
Tp = 2458153.516418316867
N  = 6.462582355082777E-01
MA = 5.556759778122834E+01
TA = 8.379865292065897E+01
A  = 1.324947455923051E+00
AD = 1.663753899326148E+00
PR = 5.570528624936787E+02
elset_time = datetime(2018, 5, 1)

# ==============================================================================

# Compute Roadster's orbit from JPL elements
Roadster = CelestialBody("Roadster")
Roadster.orbit = Orbit(primary=Sun, elements=(A*AU, EC, IN, W, OM, MA, elset_time))
Roadster_elems = Roadster.get_elements()

# Create figure
plt.style.use('dark_background')
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.set_aspect('equal')

# For English formatted dates
locale.setlocale(locale.LC_ALL, 'en_US.utf8')

# Orbits and planet positions
Mercury.plot_orbit(ax=ax, proj="2D", color=Mercury_color, units=AU)
Mercury.plot_at_time(diagram_date)
Mercury_x,Mercury_y,_,_,_,_ = Mercury.posvel_at_time(diagram_date)

Venus.plot_orbit(ax=ax, proj="2D", color=Venus_color, units=AU)
Venus.plot_at_time(diagram_date)
Venus_x,Venus_y,_,_,_,_ = Venus.posvel_at_time(diagram_date)

Earth.plot_orbit(ax=ax, proj="2D", color=Earth_color, units=AU)
Earth.plot_at_time(diagram_date)
Earth_x,Earth_y,_,_,_,_ = Earth.posvel_at_time(diagram_date)

Mars.plot_orbit(ax=ax, proj="2D", color=Mars_color, units=AU)
Mars.plot_at_time(diagram_date)
Mars_x,Mars_y,_,_,_,_ = Mars.posvel_at_time(diagram_date)

Roadster.plot_orbit(ax=ax, proj="2D", color=Roadster_color, lw=3.0, units=AU)
Roadster_x,Roadster_y,_,_,_,_ = Roadster.posvel_at_time(diagram_date)

# Sun
plt.plot([0],[0],"o", color=Sun_color)
plt.text(0, -0.17, "Sun", color=Sun_color, ha="center")

# Annotations

# Semi-transparent black backgroud to make annotations easier to read
bbox = dict(pad=0, fc=(0,0,0,0.7), ec="none", lw=0)

plt.annotate("Mercury", xy=(Mercury_x/AU, Mercury_y/AU), xytext=(-5, -23), textcoords="offset pixels", bbox=bbox, color=Mercury_color, ha="center", fontsize=12)

plt.annotate("Venus", xy=(Venus_x/AU, Venus_y/AU), xytext=(35, -6), textcoords="offset pixels", bbox=bbox, color=Venus_color, ha="center", fontsize=12)

plt.annotate("Earth", xy=(Earth_x/AU, Earth_y/AU), xytext=(35, -6), textcoords="offset pixels", bbox=bbox, color=Earth_color, ha="center", fontsize=12)

plt.annotate("Mars", xy=(Mars_x/AU, Mars_y/AU), xytext=(-30, -6), textcoords="offset pixels", bbox=bbox, color=Mars_color, ha="center", fontsize=12)

plt.annotate("Launch\n%s" % launch_date.strftime("%b %d, %Y"), xy=(Earth_x/AU, Earth_y/AU), xytext=(-65,0), textcoords="offset pixels", bbox=bbox, color=Roadster_color, ha="center", fontsize=12)

plt.annotate("Roadster\norbit", xy=(-1.08, -0.26), xytext=(-35, 20), textcoords="offset pixels", bbox=bbox, ha="center", fontsize=10, arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=Roadster_color), color=Roadster_color)

plt.annotate("Mars orbit crossing\n%s" % crossing1_date.strftime("%b %d, %Y"), xy=(-0.20, -1.49), xytext=(-10,-50), textcoords="offset pixels", ha="center", fontsize=10, arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))

plt.annotate("Mars orbit crossing\n%s" % crossing2_date.strftime("%b %d, %Y"), xy=(1.36, 0.41), xytext=(15,80), textcoords="offset pixels", ha="center", fontsize=10, arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))

plt.annotate("Aphelion %.2f AU\n%s" % (aphelion_dist, aphelion_date.strftime("%b %d, %Y")), xy=(1.23, -1.16), xytext=(20,-70), textcoords="offset pixels", ha="center", fontsize=10, arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))

plt.annotate("Perihelion %.2f AU\n%s" % (perihelion_dist, perihelion_date.strftime("%b %d, %Y")), xy=(-0.70, 0.71), xytext=(30,50), textcoords="offset pixels", ha="center", fontsize=10, arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))

plt.text(1.8, 0, "♈", ha="center", fontsize=14, color="gray")
plt.annotate("", xy=(1.93, -0.05), xytext=(-30,0), textcoords="offset pixels", arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="gray"))

s = "JPL Horizons elements (epoch: %s)\na=%.3f AU, e=%.3f, i=%.2f°, ω=%.1f°, Ω=%.1f°, M=%.2f°" % (elset_time.strftime("%Y-%B-%d"), Roadster_elems['a']/AU, Roadster_elems['e'], Roadster_elems['i'], Roadster_elems['arg'], Roadster_elems['LAN'], Roadster_elems['M0'])
plt.text(0, 1.73, s, ha="center", color="gray", fontsize=10)

# Copyright notice -- remove for Wikipedia use
plt.text(-1.95, -1.95, "Meithan West | CC-BY-SA", ha="left", fontsize=8, color="0.2")

# Plot settings
plt.xlim(-2,2)
plt.ylim(-2,2)
for tick in plt.gca().yaxis.get_major_ticks():
  tick.label.set_fontsize(9)
for tick in plt.gca().xaxis.get_major_ticks():
  tick.label.set_fontsize(9)
plt.xlabel("AU")
plt.ylabel("AU")
plt.subplots_adjust(top=0.99, bottom=0.05, left=0.110, right=0.975, hspace=0.2, wspace=0.2)

# Save figure
plt.savefig(fname)
print("Saved %s" % fname)
plt.show()

plt.close()
