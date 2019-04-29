set timestamp
set title "Bessel Log Error"
set xlabel "log10(# of intervals)"
set ylabel "log10(relative error)"
set logscale
plot "sum_out.dat" u 1:4
