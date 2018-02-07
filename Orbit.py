# A class to represent general 3D orbits, with the ability to plot them
# (using matplotlib) or compute position and/or velocity at arbitrary
# times
from datetime import datetime, timedelta
from math import pi, sin, cos, tan, atan, acos
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from CelestialBody import *
import utils

class Orbit:

  # The only required argument is the primary, which must be a
  # CelestialBody instance. The optional argument elements must be
  # contain the elements in the form specified by from_elements()
  def __init__(self, primary, elements=None):

    self.a = None
    self.e = None
    self.i = None
    self.w = None
    self.LAN = None
    self.M0 = None
    self.t0 = None
    self.long_peri = None
    self.L0 = None
    self.b = None
    self.c = None
    self.P = None
    self.n = None
    self.rpe = None
    self.rap = None
    self.vpe = None
    self.vap = None

    self.primary = primary
    if elements is not None:
      self.from_elements(*elements)

    self.ax = None
    self.proj = "3D"
    self.color = "C0"
    self.units = 1.0

  # Defines the orbit from the "standard" set of orbital elements:
  #  a: semi-major axis (meters)
  #  e: eccentricity
  #  i: inclination (degrees)
  #  arg: argument of periapsis (degrees)
  #  LAN: longitude of the ascending node (degrees)
  #  M0: mean anomaly at epoch (degrees)
  #  t0: UTC epoch of the elements, as a tz-naive datetime object
  def from_elements(self, a, e, i, arg, LAN, M0, t0):
    self.a = a
    self.e = e
    self.i = i
    self.w = arg
    self.LAN = LAN
    self.M0 = M0
    self.t0 = t0
    self.set_derived_params()

  # Defines the orbit from the 3D state vectors
  #  pos: 3D cartesian position vector, in m
  #  vel: 3D cartesian velocity vector, in m/s
  #  time: UTC time as which the state is given, tz-naive datetime object
  #  t0: (optional) UTC epoch for the elements, tz-naive datetime object
  # If t0 is not set, the orbit's epoch will be time. If t0 is given, the
  # mean anomaly at epoch will be adjusted to t0.
  def from_statevectors(self, pos, vel, time, t0=None):

    if t0 is None: t0 = time
    self.t0 = t0

    pos = np.array(pos)
    vel = np.array(vel)

    h = np.cross(pos, vel)
    ecc = np.cross(vel, h)/self.primary.mu - pos/np.linalg.norm(pos)
    n = np.cross(np.array((0,0,1)), h)
    nu = acos(np.dot(ecc,pos)/(np.linalg.norm(ecc)*np.linalg.norm(pos)))
    if np.dot(pos,vel) < 0:
      nu = 2*pi - nu
    i = acos(h[2]/np.linalg.norm(h))
    e = np.linalg.norm(ecc)
    E = 2*atan(tan(nu/2)/sqrt((1+e)/(1-e)))
    LAN = acos(n[0]/np.linalg.norm(n))
    if n[1] < 0:
      LAN = 2*pi - LAN
    w = acos(np.dot(n,ecc)/(np.linalg.norm(n)*np.linalg.norm(ecc)))
    if ecc[2] < 0:
      w = 2*pi - w
    M = E - e*sin(E)
    a = 1/(2/np.linalg.norm(pos) - np.linalg.norm(vel)**2/self.primary.mu)
    M0 = M - sqrt(self.primary.mu/a**3)*self.secs_since_epoch(time)

    i = clamp_degs(to_degs(i))
    w = clamp_degs(to_degs(w))
    LAN = clamp_degs(to_degs(LAN))
    M0 = clamp_degs(to_degs(M0))

    self.from_elements(a, e, i, w, LAN, M0, t0=t0)

  # Defines the orbit from a two-line element set
  # The argument TLE must be either a string with two newline characters
  # demarking the title and two data lines of the TLE, or a 3-tuple with
  # these elements
  def from_TLE(self, TLE):

    if isinstance(TLE, str):
      TLE = TLE.split("\n")
    assert len(TLE) == 3

    epoch_year = int(TLE[1][18:20])
    if epoch_year < 57: epoch_year += 2000
    else: epoch_year + 1900
    epoch_day = float(TLE[1][20:32])

    epoch = datetime(epoch_year, 1, 1) + timedelta(days=epoch_day-1)
    inc = float(TLE[2][8:16])
    RAAN = float(TLE[2][17:25])
    ecc = float("0."+TLE[2][26:33])
    w = float(TLE[2][34:42])
    M0 = float(TLE[2][43:51])
    motion = float(TLE[2][52:63])   # in revolutions per day

    P = 86400/motion
    a = (self.primary.mu*(P/(2*pi))**2)**(1.0/3)

    self.from_elements(a, ecc, inc, w, RAAN, M0, epoch)

  # Sets derived orbital parameters from main elements
  def set_derived_params(self):
    self.long_peri = clamp_degs(self.w + self.LAN)
    self.L0 = clamp_degs(self.M0 + self.long_peri)
    self.b = self.a*sqrt(1-self.e**2)
    self.c = self.e * self.a
    self.P = 2*pi*sqrt(self.a**3/self.primary.mu)
    self.n = 2*pi/self.P
    self.rpe = (1-self.e)*self.a
    self.rap = (1+self.e)*self.a
    self.vpe = self.vel_at_radius(self.rpe)
    self.vap = self.vel_at_radius(self.rap)

  def get_elements(self):
    return {"a": self.a, "e": self.e, "i": self.i, "arg": self.w, "LAN": self.LAN, "M0": self.M0, "long_peri": self.long_peri, "L0": self.L0, "epoch": self.t0}

  def set_axes(self, ax):
    self.ax = ax

  def set_proj(self, proj):
    self.proj = proj

  # Returns the (fractional) number of seconds since the epoch
  # Negative it before epoch
  def secs_since_epoch(self, time):
    return (time - self.t0).total_seconds()

  # Speed at given radius (i.e. vis-viva equation)
  def vel_at_radius(self, radius):
    assert self.rpe <= radius <= self.rap
    return sqrt(self.primary.mu*(2/radius - 1/self.a))

  # Radius at given true anomaly
  def radius_at_nu(self, nu):
    return self.a*(1-self.e**2)/(1+self.e*cos(utils.to_rads(nu)))

  # Mean anomaly M at provided time
  # Time must be a UTC tz-naive datetime object
  def M_at_time(self, time):
    return self.n*self.secs_since_epoch(time) + utils.to_rads(self.M0)

  # Returns the position and velocity vectors (in the orbital
  # frame) at arbitrary time (same clock as epoch)
  # Time must be a UTC tz-naive datetime object
  def posvel_at_time(self, time):

    # Get mean anomaly at specified time
    M = self.M_at_time(time)

    # Solve Kepler and obtain true anomaly
    E0 = -self.primary.mu/(2*self.a)
    if E0 < 0:
      E = self.solveKepler(M, self.e)
      nu = 2*atan(sqrt((1+self.e)/(1-self.e))*tan(E/2.0))
    else:
      pass
      # H = self.solveKepler(M, self.ECC)
      # nu = 2*atan(sqrt((self.ECC+1)/(self.ECC-1))*tanh(H/2.0))
    # if self.hz < 0: nu *= -1

    # Perifocal position
    p = self.a*(1.0-self.e**2)
    r = p/(1 + self.e*cos(nu))
    x = r*cos(nu)
    y = r*sin(nu)
    z = 0

    xx = self.a*(cos(E)-self.e)
    yy = self.a*sqrt(1-self.e**2)*sin(E)
    # print(x, xx, abs(x-xx))
    # print(y, yy, abs(y-yy))

    # Perifocal velocity
    vx = sqrt(self.primary.mu*self.a)/r*(-sin(E))
    vy = sqrt(self.primary.mu*self.a)/r*(sqrt(1-self.e**2)*cos(E))
    vz = 0

    # Transform to final reference frame
    x, y, z = self.transform(x, y, z)
    vx, vy, vz = self.transform(vx, vy, vz)

    return x, y, z, vx, vy, vz

  # Solves Kepler's equation, M = E - ecc*sin(E), for the
  # eccentric anomaly E, given M and ecc; if ecc > 0, solves the
  # hyperbolic Kepler equation instead, M = ecc*sinh(H) - H, for H.
  def solveKepler(self, M, ecc):

    # Tolerance (relative) for Newton's method
    tol = 1e-6
    max_iter = 1000

    if (M == 0):
      return 0

    while M > pi: M -= 2*pi
    while M < pi: M += 2*pi

    if ecc > 0.9:
      E = pi if M < 2*pi else -pi
    else:
      E = M

    rel_err = tol*2.0
    iter = 1
    while (rel_err > tol):
      if (ecc < 1):
        Enew = E - (E-ecc*sin(E)-M)/(1.0-ecc*cos(E))
      else:
        Enew = E - (ecc*sinh(E)-E-M)/(ecc*cosh(E)-1.0)
      if (E != 0): rel_err = abs((Enew-E)/E)
      iter += 1
      if iter > max_iter:
        print("\nsolveKepler solver reached max iterations!")
        print("M=", M)
        print("ecc=", ecc)
        print("current E=", Enew)
        print("last E=", E)
        raise Exception("Kepler solver reached max iterations!")
      E = Enew

    return E

  # Applies the coordinate transform from perifocal coordinates to
  # the final reference frame coordinates (e.g. ecliptic)
  def transform(self, x, y, z):
    # x1, y1 = rotate2D(x, y, self.w+self.LAN)
    # z1 = z
    # nx, ny = rotate2D(1, 0, self.LAN)
    # x2, y2, z2 = rotate3D(x1, y1, z1, self.i, nx, ny, 0)
    # return x2, y2, z2
    # x, y, z = rotate3D(x, y, z, -self.w, 0, 0, 1)
    # x, y, z = rotate3D(x, y, z, -self.i, 1, 0, 0)
    # x, y, z = rotate3D(x, y, z, -self.LAN, 0, 0, 1)
    x, y, z = Rz(x, y, z, self.w)
    x, y, z = Rx(x, y, z, self.i)
    x, y, z = Rz(x, y, z, self.LAN)
    return x, y, z

  # Plots the position of the body at the given time
  # Time must be a UTC tz-naive datetime object
  def plot_at_time(self, time, proj=None, ax=None, color=None, units=None):

    if ax is not None:
      self.ax = ax
    if proj is not None:
      self.proj = proj
    if color is not None:
      self.color = color
    if units is not None:
      self.units = units

    x, y, z, vx, vy, vz = self.posvel_at_time(time)
    x /= self.units; y /= self.units; z /= self.units

    if self.ax is None:
      fig = plt.figure()
      if self.proj == "3D":
        self.ax = fig.add_subplot(111, projection='3d')
      elif self.proj == "2D":
        self.ax = fig.add_subplot(111)
      self.ax.set_aspect('equal')

    if self.proj == "3D":
      self.ax.scatter([x],[y],[z], color=self.color, zorder=10)
    elif self.proj == "2D":
      self.ax.scatter([x],[y], s=10, color=self.color, zorder=10)

  # Plots the orbit
  # Can optionally provide the axes
  def plot(self, proj=None, ax=None, show_apsides=False, show_nodes=False, show_axes=False, show_primary=False, color=None, units=None):

    if ax is not None:
      self.ax = ax
    if proj is not None:
      self.proj = proj
    if color is not None:
      self.color = color
    if units is not None:
      self.units = units

    # Generate standard ellipse (perifocal coords) and transform it
    xs = []; ys = []; zs = []
    xs1, ys1 = gen_ellipse(self.a, self.b)
    xs1 = [x-self.c for x in xs1]
    for i in range(len(xs1)):
      x, y, z = self.transform(xs1[i], ys1[i], 0)
      x /= self.units; y /= self.units; z /= self.units
      xs.append(x); ys.append(y); zs.append(z)

    if self.ax is None:
      fig = plt.figure()
      if self.proj == "3D":
        self.ax = fig.add_subplot(111, projection='3d')
      elif self.proj == "2D":
        self.ax = fig.add_subplot(111)
      self.ax.set_aspect('equal')

    # Plot orbit
    if self.proj == "3D": self.ax.plot(xs, ys, zs, color=self.color, lw=1.0)
    elif self.proj == "2D": self.ax.plot(xs, ys, color=self.color, lw=1.0)

    scale = self.rap*1.2 / self.units
    self.ax.set_xlim(-scale, scale)
    self.ax.set_ylim(-scale, scale)
    if self.proj == "3D":
      self.ax.set_zlim(-scale, scale)

    if show_primary:

      # Primary as a circle/sphere
      R = self.primary.radius/self.units
      if self.proj == "3D":
        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 15)
        x = R * np.outer(np.cos(u), np.sin(v))
        y = R * np.outer(np.sin(u), np.sin(v))
        z = R * np.outer(np.ones(np.size(u)), np.cos(v))
        self.ax.plot_surface(x, y, z, rstride=1, cstride=1, color="gray")
      elif self.proj == "2D":
        import matplotlib.patches as mpatches
        circle = mpatches.Circle((0,0), R, fill=True, color="gray")

    # else:
    #
    #   # Just at mark at the center of gravity
    #   if self.proj == "3D":
    #     self.ax.plot([0], [0], [0], "x", color="gray")
    #   elif self.proj == "2D":
    #     self.ax.plot([0], [0], "x", color="gray")

    if show_apsides:

      # Position and line of apsides
      xpe, ype, zpe = self.transform(self.rpe, 0, 0)
      xap, yap, zap = self.transform(-self.rap, 0, 0)
      xpe /= self.units; ype /= self.units; zpe /= self.units
      xap /= self.units; yap /= self.units; zap /= self.units
      if self.proj == "3D":
        self.ax.plot([xpe], [ype], [zpe], "ko")
        self.ax.plot([xap], [yap], [zap], "ko")
        self.ax.plot([xpe,xap], [ype,yap], [zpe,zap], "k", ls=(0,(3,3)))
      elif self.proj == "2D":
        self.ax.plot([xpe], [ype], "ko")
        self.ax.plot([xap], [yap], "ko")
        self.ax.plot([xpe,xap], [ype,yap], "k", ls=(0,(3,3)))

    if show_nodes:

      # Nodes and line of nodes
      r1 = self.radius_at_nu(-self.w)
      xno1 = r1
      yno1 = 0
      xno1, yno1 = rotate2D(xno1, yno1, self.LAN)
      r2 = self.radius_at_nu(180-self.w)
      xno2 = -r2
      yno2 = 0
      xno2, yno2 = rotate2D(xno2, yno2, self.LAN)
      xno1 /= self.units; yno1 /= self.units
      xno2 /= self.units; yno2 /= self.units
      self.ax.plot([xno1], [yno1], "bo")
      self.ax.plot([xno2], [yno2], "bo")
      self.ax.plot([xno1,xno2], [yno1,yno2], "b", ls=(3,(3,3)))

    if show_axes:

      # Coordinate axes
      self.ax.plot([0,self.a/self.units], [0,0], "r")
      self.ax.plot([0,0], [0,self.a/self.units], "g")
      if self.proj == "3D":
        self.ax.plot([0,0], [0,0], [0,self.a/self.units], "b")
