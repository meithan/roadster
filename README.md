# roadster
Code to make an animation of the orbit of the Roadster to be launched by the Falcon Heavy:

https://www.youtube.com/watch?v=i8XDwt-67vY

Requires the Python plotting library [matplotlib](https://matplotlib.org/users/installing.html) (and numpy, which should come with matplotlib).

To make the animation:

1) Run Roadster.py to generate the animation's frames (the images will be saved in the 'frames' folder)
2) In Linux, one can use ffmpeg to assemble a video from the frames (see ffmpeg.txt for the command)

Should work in both Python 2.x and Python 3.


