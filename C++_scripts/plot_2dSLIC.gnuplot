set terminal gif animate delay 5 size 1000, 800 background rgb "#d3d3d3"
set output 'lightning_density.gif'


set view 43, 321  # diagonal from above

# Remove axes and labels


# lock range for colorbar across frames

set xrange [-1:1]
set yrange [-1:1]
set zrange [0:10]

set title "Density Evolution for Cylindrical Sod Test"

do for [i=0:79] {
    splot sprintf('SLIC%d.dat', i) using 1:2:3 
}
