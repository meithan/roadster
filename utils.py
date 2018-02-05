# Generic utility functions and constants

from math import pi, sqrt, sin, cos
from datetime import datetime

# ==================================

# Standard gravitational acceleration (m/s^2)
g0 = 9.80665

# 2014 CODATA-recommended value (m^3 kg^-1 s^-2)
GRAV = 6.67408e-11

# Astronomical constants
AU = 149597870700.0

# J2000 epoch
J2000 = datetime(2000,1,1,11,58,56)

# ==================================

# Converts angle in degrees to radians
def to_rads(degrees):
  return degrees / 180. * pi

# Converts andle in radians to degrees:
def to_degs(radians):
  return radians / pi * 180.

def clamp_degs(degrees):
  while degrees > 360:
    degrees -= 360
  while degrees < 0:
    degrees += 360
  return degrees

# Generates a "standard" planar ellipse with the given semi-major
# and semi-minor axes
def gen_ellipse(a, b):

  xs = []
  ys = []
  step = 0.01
  theta = 0
  while theta < 2*pi:
    xs.append(a*cos(theta))
    ys.append(b*sin(theta))
    theta += step

  return xs, ys

# Planar rotation of the given point through angle theta (degrees)
def rotate2D(x, y, theta):
  theta_rad = to_rads(theta)
  xp = x*cos(theta_rad) - y*sin(theta_rad)
  yp = x*sin(theta_rad) + y*cos(theta_rad)
  return xp, yp

# 3D rotation of the given point through angle theta (degrees) around
# axis specified by ux, uy, uz
def rotate3D(x, y, z, theta, ux, uy, uz):
  if (ux**2 + uy**2 + uz**2 - 1) > 1e-3:
    norm = sqrt(ux**2 + uy**2 + uz**2)
    ux /= norm; uy /= norm; uz /= norm
  theta_rad = to_rads(theta)
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = x*(C+ux**2*(1-C)) + y*(ux*uy*(1-C)-uz*S) + z*(ux*uz*(1-C)+uy*S)
  yp = x*(uy*ux*(1-C)+uz*S) + y*(C+uy**2*(1-C)) + z*(uy*uz*(1-C)-ux*S)
  zp = x*(uz*ux*(1-C)-uy*S) + y*(uz*uy*(1-C)+ux*S) + z*(C+uz**2*(1-C))
  return xp, yp, zp

# 3D elementary rotation around the x-axis
def Rx(x, y, z, theta):
  theta_rad = to_rads(theta)
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = x
  yp = C*y - S*z
  zp = S*y + C*z
  return xp, yp, zp

# 3D elementary rotation around the y-axis
def Ry(x, y, z, theta):
  theta_rad = to_rads(theta)
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = C*x + S*z
  yp = y
  zp = -S*x + C*z
  return xp, yp, zp

# 3D elementary rotation around the z-axis
def Rz(x, y, z, theta):
  theta_rad = to_rads(theta)
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = C*x - S*y
  yp = S*x + C*y
  zp = z
  return xp, yp, zp
