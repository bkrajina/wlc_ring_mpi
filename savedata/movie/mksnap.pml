delete *
reset

cd pdb

movie.load snap*.pdb, snap


hide all

set stick_radius=1.0
show spheres, (name A1)
color blue, (name A1)
alter (name A1), vdw=0.05
show sticks, (name A1)
rebuild

#run ../color_b.py
#color_b snap,mode=hist,gradient=bgr,minimum=0,maximum=1,nbins=40,sat=1.,value=1.

bg_color white
reset

cd ..

#cd ../png
#set ray_trace_frames = 1
#mpng snap