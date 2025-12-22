# GNUPLOT script for trajectory comparison
set terminal pngcairo size 1400,900 enhanced font 'Arial,14'
set output 'IL76_comparison.png'

set title 'Ил-76: Сравнение оптимальных траекторий'
set xlabel 'Скорость, км/ч'
set ylabel 'Высота, м'
set grid
set key top left box
set datafile separator ','

plot 'il76_time_optimal_full.csv' every ::1 using 5:3 \
     with lines lw 3 lc rgb 'red' \
     title 'Минимум времени', \
     'il76_fuel_optimal_full.csv' every ::1 using 5:3 \
     with lines lw 3 lc rgb 'blue' \
     title 'Минимум топлива', \
     'il76_time_optimal_full.csv' every 50::1 using 5:3 \
     with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     'il76_fuel_optimal_full.csv' every 50::1 using 5:3 \
     with points pt 7 ps 1.5 lc rgb 'blue' notitle
