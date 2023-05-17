ef=  6.00000000
set xtics ( \
"M"   -0.334813, \
"{/Symbol \107}"    0.000000, \
"K"    0.289957 )
set xrange [   -0.334813:    0.289957]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 7.00in, 7.50in
set output "band_checks.pdf"
unset colorbox
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set ylabel "E-E_{CBM} (eV)"
set yrange [ -1 : 1.0 ]
unset key
set ytics 1.0 scale 1 nomirror out
set mytics 2
set parametric
set trange [-10:10]
plot "band.dat" u 1:2:3 with image