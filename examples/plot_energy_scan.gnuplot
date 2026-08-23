# Plot the p-11B S-factor and cross section in the style of Figure 1.
#
# Use the benchmark directly:
#   gnuplot -c examples/plot_energy_scan.gnuplot
#
# Or plot an already saved scan and optionally choose the output filename:
#   ./build/benchmark_pb11 --energy-scan > energy_scan.csv
#   gnuplot -c examples/plot_energy_scan.gnuplot energy_scan.csv pb11_figure1.png

input_source = "< ./build/benchmark_pb11 --energy-scan"
output_file = "pb11_figure1.png"

if (ARGC >= 1 && strlen(ARG1) > 0) input_source = ARG1
if (ARGC >= 2 && strlen(ARG2) > 0) output_file = ARG2

set terminal pngcairo enhanced size 1200,1000 font "DejaVu Sans,18"
set output output_file

set datafile separator comma
set border linewidth 1.5
set tics out nomirror
set grid xtics ytics mxtics mytics back lc rgb "#d0d0d0" lw 0.8
set key off

set logscale x 10
set xrange [0.1:10]
set mxtics 10

# Published piece boundaries, shown in both panels.
set arrow 1 from 0.400, graph 0 to 0.400, graph 1 nohead \
    dt 2 lw 1.2 lc rgb "#777777"
set arrow 2 from 0.668, graph 0 to 0.668, graph 1 nohead \
    dt 2 lw 1.2 lc rgb "#777777"

set multiplot layout 2,1 rowsfirst \
    margins 0.12,0.96,0.10,0.96 spacing 0.0,0.07

# Panel (a): logarithmic S-factor axis resolves the narrow 148 keV resonance
# and the several-orders-of-magnitude high-energy decrease.
set logscale y 10
set yrange [0.1:5000]
set mytics 10
set format x ""
set format y "10^{%L}"
set ylabel "S(E)  [MeV barn]" offset 1.0,0
set label 1 "(a)" at graph 0.02,0.92 front
plot input_source using 1:2 with lines lw 2.5 lc rgb "#e41a1c"

# Panel (b): the cross section is linear in y, matching Figure 1(b).
unset logscale y
set yrange [0:1.5]
set ytics 0.2
set mytics 2
set xtics ("10^{-1}" 0.1, "10^{0}" 1.0, "10^{1}" 10.0)
set format x "%g"
set format y "%.1f"
set xlabel "Centre-of-mass energy E  [MeV]"
set ylabel "{/Symbol s}(E)  [barn]" offset 1.0,0
set label 1 "(b)" at graph 0.02,0.92 front
plot input_source using 1:3 with lines lw 2.5 lc rgb "#e41a1c"

unset multiplot
unset output

print sprintf("Wrote %s", output_file)
