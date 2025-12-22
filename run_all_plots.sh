#!/bin/bash
echo "Creating IL-76 trajectory plots..."
echo
echo "1. Comparison plot (V-H diagram)"
gnuplot plot_comparison.gp
echo "2. Altitude vs time"
gnuplot plot_altitude_time.gp
echo "3. Speed vs time"
gnuplot plot_speed_time.gp
echo "4. Fuel consumption"
gnuplot plot_fuel_time.gp
echo
echo "All plots created:"
echo "- IL76_comparison.png"
echo "- IL76_altitude_time.png"
echo "- IL76_speed_time.png"
echo "- IL76_fuel_time.png"
chmod +x run_all_plots.sh
