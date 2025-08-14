set terminal pngcairo size 1000,600
set output "lightning_selected.png"


plot \
  "wrapped_0.dat" using 1:2 title "t = 0μs" with lines, \
  "wrapped_1.dat" using 1:2 title "t = 10μs" with lines, \
  "wrapped_2.dat" using 1:2 title "t = 20μs" with lines, \
  "wrapped_4.dat" using 1:2 title "t = 40μs" with lines, \
  "wrapped_6.dat" using 1:2 title "t = 60μs" with lines, \
  "wrapped_10.dat" using 1:2 title "t = 100μs" with lines, \
  "wrapped_15.dat" using 1:2 title "t = 150μs" with lines
