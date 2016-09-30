#!/bin/env ruby

input_file=ARGV[0]
output_file=input_file.sub(/\..*$/,'.png')

# Get model forcing and delta_t
amplitude=0.0
bias=0.0
frequency=0.0
phase=0.0
delta_t=0.0
f = File.new(input_file)
f.each do |line|
  break if line=~/^\s*$/
  if line=~/^model_amplitude,\s*(\S+)$/
    amplitude=$1
  end
  if line=~/^model_bias,\s*(\S+)$/
    bias=$1
  end
  if line=~/^model_frequency,\s*(\S+)$/
    frequency=$1
  end
  if line=~/^model_phase,\s*(\S+)$/
    phase=$1
  end
  if line=~/^model_delta_t,\s*(\S+)$/
    delta_t=$1
  end
end
f.close

# Create a plot for each reservation
IO.popen("gnuplot","w") do |gnuplot|
  gnuplot.puts "set datafile separator ','"
  gnuplot.puts "set xlabel('X')"
  gnuplot.puts "set ylabel('State(x)')"
  gnuplot.puts "set title 'Sine Model State (amplitude=#{amplitude}, bias=#{bias}, frequency=#{frequency}, phase=#{phase}, delta_t=#{delta_t})'"
  gnuplot.puts "set output '#{output_file}'"
  gnuplot.puts "set terminal png size 1000,600"
  gnuplot.puts "plot '#{input_file}' index 1 using 1:3 with linespoints"
end

