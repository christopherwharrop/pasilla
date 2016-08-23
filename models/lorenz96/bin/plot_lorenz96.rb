#!/bin/env ruby

input_file=ARGV[0]
output_file=input_file.sub(/\..*$/,'.png')

# Create a plot for each reservation
IO.popen("gnuplot","w") do |gnuplot|
  gnuplot.puts "set datafile separator ','"
  gnuplot.puts "set xlabel('X')"
  gnuplot.puts "set ylabel('State(x)')"
  gnuplot.puts "set title 'Lorenz96 Model State'"
  gnuplot.puts "set output '#{output_file}'"
  gnuplot.puts "set terminal png size 1000,600"
  gnuplot.puts "plot '#{input_file}' index 1 using 1:3 with linespoints"
end

