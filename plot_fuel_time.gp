# Fuel consumption
set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'IL76_fuel_time.png'

set title 'Ил-76: Расход топлива (ДП)'
set xlabel 'Время, с'
set ylabel 'Израсходованное топливо, кг'
set grid
set key top left
set datafile separator ','

plot 'il76_dp_solution.csv' using 1:9ё with lines lw 2 lc rgb 'purple' title 'Топливо'
