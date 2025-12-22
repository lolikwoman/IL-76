@echo off
echo Построение графиков Ил-76 (ДП)...
echo.
gnuplot plot_VH.gp
gnuplot plot_altitude_time.gp
gnuplot plot_speed_time.gp
gnuplot plot_fuel_time.gp
echo.
echo Графики построены:
echo  - IL76_VH.png
echo  - IL76_altitude_time.png
echo  - IL76_speed_time.png
echo  - IL76_fuel_time.png
pause
