from utils import *

# A celestial body
class CelestialBody:

  def __init__(self, name, mu=None, radius=None, rot_period=None, orbit=None):
    self.name = name
    self.mu = mu
    if mu is not None:
      self.mass = GRAV/self.mu
    self.radius = radius
    self.rot_period = rot_period
    self.orbit = orbit
    self.color = 'C0'

  def get_elements(self):
    return self.orbit.get_elements()

  def plot_at_time(self, time, ax=None, color=None, proj=None, units=None):
    if color is not None:
      self.color = color
    return self.orbit.plot_at_time(time, ax=ax, color=self.color, proj=proj, units=units)

  def plot_orbit(self, proj=None, ax=None, show_apsides=False, show_nodes=False, show_axes=False, show_primary=False, color=None, units=None, lw=None):
    if color is not None:
      self.color = color
    self.orbit.plot(proj=proj, ax=ax, show_apsides=show_apsides, show_nodes=show_nodes, show_axes=show_axes, show_primary=show_primary, color=self.color, units=units, lw=lw)

  def posvel_at_time(self, time):
    return self.orbit.posvel_at_time(time)
