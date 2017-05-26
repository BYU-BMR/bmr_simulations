import ovito
from ovito.vis import Viewport,RenderSettings,TachyonRenderer
from ovito.io import import_file
import sys,re,math

filename = "moving_hopper.lammpstrj"

for arg in sys.argv:
	print(arg)
	if re.match(".*.lammpstrj",arg):
		filename = arg



print("Opening %s using OVITO %i.%i.%i" % ((filename,) + ovito.version))


node = import_file(filename, multiple_frames = True)
print("Num Frames:",node.source.num_frames)

node.add_to_scene()
cell = node.source.cell
centerPoint = (
	(cell.vector1[0]-cell.origin[0])/2.0,
	(cell.vector2[1]-cell.origin[1])/2.0,
	(cell.vector3[2]-cell.origin[2])/2.0,
	)
cameraPosition = (centerPoint[0],-800,centerPoint[2])

#frameToLoad = node.source.num_frames - 1
frameToLoad = math.floor(node.source.num_frames*0.97)
ovito.dataset.anim.current_frame = frameToLoad
print("Loaded frame", ovito.dataset.anim.current_frame+1)

image_rs = RenderSettings(
    filename = 'last_frame.png',
    size = (1024,768),
    background_color = (0.8,0.8,1.0),
    renderer = TachyonRenderer()
)

movie_rs = RenderSettings(
	filename = 'movie.avi',
    size = (1024,768),
    background_color = (0.8,0.8,1.0),
    renderer = TachyonRenderer(),
    #range = RenderSettings.Range.ANIMATION,
    range = RenderSettings.Range.CUSTOM_INTERVAL,
    custom_range = (0,850)
)

vp = Viewport()
vp.type = Viewport.Type.PERSPECTIVE
vp.camera_pos = cameraPosition
vp.camera_dir = (0, 1, 0)

vp.render(movie_rs)

#ovito.dataset.viewports.active_vp.render(rs)
