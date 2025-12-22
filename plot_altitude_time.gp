# Altitude vs time
set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'IL76_altitude_time.png'

set title 'Ил-76: Высота от времени (ДП)'
set xlabel 'Время, с'
set ylabel 'Высота, м'
set grid
set key top left
set datafile separator ','

plot 'il76_dp_solution.csv' using 1:3 with lines lw 2 lc rgb 'dark-green' title 'Высота'
