using Dates

stime=DateTime(0,7,1)

# time
ΔT=30; #seconds
ndays=100; # duration of simulation in days

n=Int(round(ndays*24*60*60/ΔT)); # number of time steps

