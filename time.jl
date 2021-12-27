using Dates

stime=DateTime(0,7,1)

# time
ΔT=30; #seconds
ndays=2; # duration of simulation in days (e.g. 2 for testing/dev; 100 otherwise)

n=Int(round(ndays*24*60*60/ΔT)); # number of time steps

