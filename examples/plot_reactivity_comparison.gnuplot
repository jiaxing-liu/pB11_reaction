# Compare the direct Maxwellian integral with the independent analytic fit.
#
# Run directly from the benchmark:
#   gnuplot -c examples/plot_reactivity_comparison.gnuplot
#
# Or use an already saved scan and choose the output filename:
#   ./build/benchmark_pb11 --reactivity-scan > reactivity_scan.csv
#   gnuplot -c examples/plot_reactivity_comparison.gnuplot \
#       reactivity_scan.csv pb11_reactivity_comparison.png

input_source = "< ./build/benchmark_pb11 --reactivity-scan"
output_file = "pb11_reactivity_comparison.png"

if (ARGC >= 1 && strlen(ARG1) > 0) input_source = ARG1
if (ARGC >= 2 && strlen(ARG2) > 0) output_file = ARG2

set terminal pngcairo enhanced size 1200,1000 font "DejaVu Sans,18"
set output output_file

set datafile separator comma
set border linewidth 1.5
set tics out nomirror
set grid xtics ytics mxtics mytics back lc rgb "#d0d0d0" lw 0.8

set logscale x 10
set xrange [10:500]
set mxtics 10

# LT/HT analytic-fit boundary.
set arrow 1 from 70, graph 0 to 70, graph 1 nohead \
    dt 3 lw 1.4 lc rgb "#666666"

set multiplot layout 2,1 rowsfirst \
    margins 0.13,0.96,0.10,0.96 spacing 0.0,0.07

# Top panel: the two independently calculated reactivities.
set logscale y 10
set yrange [1e-28:1e-21]
set mytics 10
set format x ""
set format y "10^{%L}"
set ylabel "<{/Symbol s}v>  [m^{3} s^{-1}]" offset 1.2,0
set key at graph 0.97,0.08 right bottom opaque box
set label 1 "(a)" at graph 0.02,0.92 front
plot input_source using 1:2 with lines lw 3.0 lc rgb "#1f3b73" \
         title "Eq. (6) numerical integral", \
     input_source using 1:3 with lines dt 2 lw 3.0 lc rgb "#e41a1c" \
         title "Eq. (7)-(15) analytic fit"

# Bottom panel: signed relative error makes the agreement and 70 keV switch
# visible even where the two curves overlap in the logarithmic top panel.
unset logscale y
set yrange [-2.5:2.5]
set ytics 1.0
set mytics 2
set format x "10^{%L}"
set format y "%.1f"
set xlabel "Ion temperature kT  [keV]"
set ylabel "Relative error  [%]" offset 0.5,0
unset key
set label 1 "(b)" at graph 0.02,0.90 front
set arrow 2 from graph 0, first 0 to graph 1, first 0 nohead \
    lw 1.2 lc rgb "#333333"
plot input_source using 1:($4*100.0) with linespoints \
    pt 7 ps 0.35 lw 2.0 lc rgb "#2b8c4b"

unset multiplot
unset output

print sprintf("Wrote %s", output_file)
