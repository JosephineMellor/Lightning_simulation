set terminal gif size 1000, 800 animate delay 50

set output 'O_T.gif'

set pm3d
unset surface
set view map 
set cbrange[0:7]
set xrange[0:1]
set yrange[0:1]
do for [idx =0:6]{\
splot 'MHD'.idx.'.dat' u 1:2:3 w l notitle}