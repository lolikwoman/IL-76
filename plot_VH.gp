# V-H diagram
set terminal pngcairo size 1400,900 enhanced font 'Arial,14'
set output 'IL76_VH.png'

set title 'Ил-76: Диаграмма скорость–высота (ДП)'
set xlabel 'Скорость, км/ч'
set ylabel 'Высота, м'
set grid
set key top left
set datafile separator ','

plot 'il76_dp_solution.csv' using 3:2 with lines lw 3 lc rgb 'blue' title 'Оптимальная траектория'
