set timestamp
set title "Gorund State Harmonic Oscillator"
set xlabel "r"
set ylabel "Psi(r)"
plot [0:10] 'eigen_basis_dim=1_b=2.5.dat' u 1:2, 'eigen_basis_dim=5_b=2.5.dat' u 1:2,'eigen_basis_dim=10_b=2.5.dat' u 1:2,'eigen_basis_dim=20_b=2.5.dat' u 1:2, 2*x*exp(-x)
