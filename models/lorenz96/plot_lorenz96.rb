#!/bin/env ruby

# Create a plot for each reservation
IO.popen("gnuplot","w") do |gnuplot|
  gnuplot.puts "set xlabel('t')"
  gnuplot.puts "set ylabel('x(t)')"
  gnuplot.puts "set title 'Lorenz96 Trajectory'"
  gnuplot.puts "set output 'lorenz96.png'"
  gnuplot.puts "set terminal png size 1000,600"
  gnuplot.puts "plot 'lorenz.out' with linespoints"
end

