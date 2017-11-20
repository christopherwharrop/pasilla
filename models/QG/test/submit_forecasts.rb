#!/usr/bin/ruby

sigma_start=ARGV[0].to_f
sigma_end=ARGV[1].to_f
increment=ARGV[2].to_f
#parm1_start=ARGV[0].to_f
#parm1_end=ARGV[1].to_f
#parm1_increment=ARGV[2].to_f
#parm2_start=ARGV[3].to_f
#parm2_end=ARGV[4].to_f
#parm2_increment=ARGV[5].to_f

sigma_start.step(sigma_end, increment) do |sigma|
#parm1_start.step(parm1_end, parm1_increment) do |parm1|
#parm2_start.step(parm2_end, parm2_increment) do |parm2|
  sigmastr=sprintf("%6.4f",sigma)
#  parm1str=sprintf("%6.5f",parm1)
#  parm2str=sprintf("%6.5f",parm2)
#  system("qsub -A gsd-hpcs -lprocs=24,walltime=2:00:00 -F '0 1440 #{parm1str} #{parm2str}' ./run_qg_forecasts.sh")
  system("qsub -A gsd-hpcs -lprocs=24,walltime=2:00:00 -N run_qg_forecasts_#{sigmastr} -j oe -o run_qg_forecasts_#{sigmastr} -F '0 1440 #{sigmastr}' ./run_qg_forecasts.sh")
#end
end

