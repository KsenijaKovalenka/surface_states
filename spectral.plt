ef=  6.00000000
set xtics ( \
"M"   -0.324813, \
"{/Symbol \107}"    -0.0200000, \
"K"    0.279957 )
set xrange [   -0.334813:    0.289957]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 7.00in, 7.50in
set output "spectral.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
unset colorbox
set ylabel "E-E_{f} (eV)"
set yrange [ -1 : 1.0 ]
unset key
set ytics 0.5 scale 1 nomirror out
set mytics 2
set parametric
set trange [-10:10]
set arrow from -0.05,-1.07 to -0.3,-1.07 lw 3
set arrow from 0.01,-1.07 to 0.25,-1.07 lw 3
set palette rgbformulae 30,31,32
plot "band.dat" u 1:2:4 with image,\
