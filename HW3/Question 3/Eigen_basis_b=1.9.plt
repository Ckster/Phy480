set timestamp
set title "Gorund State Harmonic Oscillator"
set xlabel "r"
set ylabel "Psi(r)"
plot 'eigen_basis_dim=1_b=1.9.dat' u 1:2, 'eigen_basis_dim=5_b=1.9.dat' u 1:2,'eigen_basis_dim=10_b=1.9.dat' u 1:2,'eigen_basis_dim=20_b=1.9.dat' u 1:2, 2*x*exp(-x)
