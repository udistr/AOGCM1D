using Dates

stime=DateTime(0,7,1)

# time
ΔT=10; #seconds
ndays=10; # duration of simulation in days

n=Int(round(ndays*24*60*60/ΔT)); # number of time steps

