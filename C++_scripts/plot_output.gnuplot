set terminal pngcairo size 1000,600
set output "lightning_with_radiation.png"

plot "with_thermal.dat" using 1:2 with lines