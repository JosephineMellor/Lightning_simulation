set terminal gif animate delay 5 size 1000, 800 background rgb "#d3d3d3"
set output 'euler_density.gif'

set pm3d
set hidden3d
unset surface
set view 43, 321  # diagonal from above

# Remove axes and labels
unset border
unset xtics
unset ytics
unset ztics
unset xlabel
unset ylabel
unset zlabel
set palette rgbformulae 21,22,23

# lock range for colorbar across frames
set cbrange [0:1.0]
set xrange [-1:1]
set yrange [-1:1]
set zrange [0:1]

set title "Density Evolution for Cylindrical Sod Test"

do for [i=0:79] {
    splot sprintf('wrapped_%d.dat', i) using 1:2:3 with pm3d notitle
}

