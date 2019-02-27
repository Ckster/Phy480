set timestamp
set title "Double-Precision Milne Integration of sqrt(1+x^2)"
set xlabel "log10(# of intervals)"
set ylabel "log10(relative error)"
f(x)=m*x+b
fit [0:1.3] f(x) 'integ_hw2_Milne.dat' u 1:2 via m,b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',m, b)
plot "integ_hw2_Milne.dat" u 1:2, f(x) t title_f(m,b)
