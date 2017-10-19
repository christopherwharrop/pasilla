#!/usr/bin/ruby

# Get the method to optimize from argument
method=ARGV[0].to_i

# Get the resolution of sigma
dsigma=ARGV[1].to_f

# Get the resolution of alpha
dalpha=ARGV[2].to_f

# Set bounds of sigma
lsigma=1.0
usigma=4.0

# Set bounds of alpha
lalpha=dalpha
ualpha=0.1

# Detect number of cores on the machine
ncores=`cat /proc/cpuinfo | grep processor | wc -l`.to_i

# Set up parameter space
params=[]
lsigma.step(usigma, dsigma) do |sigma|
  lalpha.step(ualpha, dalpha) do |alpha|
    params.push({"sigma" => sigma, "alpha" => alpha})
  end
end

while !params.empty?
  batch=params.shift(ncores)
  threads=[]
  batch.each do |param|
    printf "Processing method = %d, sigma = %5.2f, alpha=%7.4f\n", method, param["sigma"], param["alpha"]
    threads << Thread.new { `./run_varsolver.sh #{method} #{param["sigma"]} #{param["alpha"]}` }
#break
  end
  threads.each do |thread|
    thread.join
  end
#break
end
