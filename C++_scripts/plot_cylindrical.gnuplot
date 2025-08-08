set terminal pngcairo size 800,600
set output 'density_surface.png'

set xlabel "x = r cos(θ)"
set ylabel "y = r sin(θ)"
set zlabel "Density"
set cbrange[0,1]
set zrange[0,1]
set hidden3d
set pm3d
set view 43, 321  # adjust for best 3D angle

splot 'wrapped.dat' with pm3d notitle

