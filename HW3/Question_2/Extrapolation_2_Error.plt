set timestamp
set title "Log Plot of 2nd Richardson Extrapolation Error"
set xlabel "log10(h)"
set ylabel "log10(relative error)"
f(x)=m*x+b
fit [-1:0] f(x) 'derivative_test.dat' u 1:5 via m,b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',m, b)
plot 'derivative_test.dat' u 1:5, f(x) t title_f(m,b)
