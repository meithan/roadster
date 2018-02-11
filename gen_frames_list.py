# Generate a file list for ffmpeg with "frozen" frames

def gen_fname(frame):
  return "frame%04i.png" % frame

min_frame = 0
max_frame = 1700
freeze_frames = [0, 156, 275, 415, 555, 974, 1147]
freeze_length = 60

frame = min_frame
while frame <= max_frame:
  if frame in freeze_frames:
    for i in range(freeze_length):
      print("file '%s'" % gen_fname(frame))
  else:
    print("file '%s'" % gen_fname(frame))
  frame += 1
