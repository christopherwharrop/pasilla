#!/usr/bin/ruby

#sigma_start=ARGV[0].to_f
#sigma_end=ARGV[1].to_f
#increment=ARGV[2].to_f
parm1_start=ARGV[0].to_f
parm1_end=ARGV[1].to_f
parm1_increment=ARGV[2].to_f
parm2_start=ARGV[3].to_f
parm2_end=ARGV[4].to_f
parm2_increment=ARGV[5].to_f

#sigma_start.step(sigma_end, increment) do |sigma|
parm1_start.step(parm1_end, parm1_increment) do |parm1|
parm2_start.step(parm2_end, parm2_increment) do |parm2|
#  system("qsub -A gsd-hpcs -lprocs=24,walltime=00:10:00 -F '#{sigmastr}' ./make_B_qg.sh")
#  sigmastr=sprintf("%6.4f",sigma)
  parm1str=sprintf("%6.5f",parm1)
  parm2str=sprintf("%6.5f",parm2)
  system("qsub -A gsd-hpcs -lprocs=24,walltime=00:10:00 -F '#{parm1str} #{parm2str}' ./make_B_qg.sh")
end
end

