# roadster
Code to make an annotated plot or an animation of the orbit of the Roadster that was launched by the Falcon Heavy:

Annotated plot: https://github.com/meithan/roadster/blob/master/Roadster_orbit.png

Animation: https://www.youtube.com/watch?v=phwywxuAbqM

Requires the Python plotting library [matplotlib](https://matplotlib.org/users/installing.html) (and numpy, which should come with matplotlib). The orbital data is included in the project.

To make an annotated static plot of the orbit:

1) Run Roadster_orbit.py (the image Roadster_orbit.png will be saved)

To produce the frames to make an animation:

1) Run Roadster_animation.py (the images will be saved in the 'frames' folder)
2) In Linux, one can use ffmpeg to assemble a video from the frames (see ffmpeg.txt for the command)

Should work in both Python 2.x and Python 3.


