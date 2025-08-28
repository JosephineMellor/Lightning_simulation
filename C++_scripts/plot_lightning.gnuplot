set terminal pngcairo size 600,600
set output "density.png"
set title "Density [kg/m^3]"


plot \
  "wrapped_1.dat"  using 1:2 title "t = 0μs"     with lines, \
  "wrapped_2.dat"  using 1:2 title "t = 10μs"    with lines, \
  "wrapped_3.dat"  using 1:2 title "t = 20μs"    with lines, \
  "wrapped_5.dat"  using 1:2 title "t = 40μs"    with lines, \
  "wrapped_7.dat"  using 1:2 title "t = 60μs"    with lines, \
  "wrapped_11.dat" using 1:2 title "t = 100μs"   with lines, \
  "wrapped_16.dat" using 1:2 title "t = 150μs"   with lines


set terminal pngcairo size 600,600
set output "temperature.png"
set title "Temperature [Kelvin]"

plot \
  "wrapped_1.dat"  using 1:5 title "t = 0μs"     with lines, \
  "wrapped_2.dat"  using 1:5 title "t = 10μs"    with lines, \
  "wrapped_3.dat"  using 1:5 title "t = 20μs"    with lines, \
  "wrapped_5.dat"  using 1:5 title "t = 40μs"    with lines, \
  "wrapped_7.dat"  using 1:5 title "t = 60μs"    with lines, \
  "wrapped_11.dat" using 1:5 title "t = 100μs"   with lines, \
  "wrapped_16.dat" using 1:5 title "t = 150μs"   with lines


set terminal pngcairo size 600,600
set output "pressure.png"
set title "Pressure [atm]"


plot \
  "wrapped_1.dat"  using 1:($4/101325) title "t = 0μs"     with lines, \
  "wrapped_2.dat"  using 1:($4/101325) title "t = 10μs"    with lines, \
  "wrapped_3.dat"  using 1:($4/101325) title "t = 20μs"    with lines, \
  "wrapped_5.dat"  using 1:($4/101325) title "t = 40μs"    with lines, \
  "wrapped_7.dat"  using 1:($4/101325) title "t = 60μs"    with lines, \
  "wrapped_11.dat" using 1:($4/101325) title "t = 100μs"   with lines, \
  "wrapped_16.dat" using 1:($4/101325) title "t = 150μs"   with lines


set terminal pngcairo size 600,600
set output "pressure2.png"


plot \
  "wrapped_1.dat"  using 1:4 title "t = 0μs"     with lines, \
  "wrapped_2.dat"  using 1:4 title "t = 10μs"    with lines, \
  "wrapped_3.dat"  using 1:4 title "t = 20μs"    with lines, \
  "wrapped_5.dat"  using 1:4 title "t = 40μs"    with lines, \
  "wrapped_7.dat"  using 1:4 title "t = 60μs"    with lines, \
  "wrapped_11.dat" using 1:4 title "t = 100μs"   with lines, \
  "wrapped_16.dat" using 1:4 title "t = 150μs"   with lines