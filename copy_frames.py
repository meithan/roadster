# Copies frames to add "frozen" frames
from shutil import copy

def gen_fname(frame):
  return "frame%04i.png" % frame

min_frame = 0
max_frame = 1700
freeze_frames = [0, 156, 276, 415, 555, 974, 1147]
freeze_length = 90
outdir = "tmp/"

# ==================

if not outdir.endswith("/"): outdir += "/"

in_frame = min_frame
out_frame = in_frame
while in_frame <= max_frame:
  in_fname = gen_fname(in_frame)
  repeat = freeze_length if in_frame in freeze_frames else 1
  for i in range(repeat):
    out_path = outdir + gen_fname(out_frame)
    copy(in_fname, out_path)
    print(out_path)
    out_frame += 1

  in_frame += 1
