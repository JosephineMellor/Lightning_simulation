set output "lightning_data.png"

do for [i=0:14] {
    plot sprintf('wrapped_%d.dat', i) using 1:2 notitle,
}

plot sprintf('wrapped_15.dat', i) using 1:2 notitle