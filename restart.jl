include("delta.jl")
include("heaviside.jl")
include("angle_of_incidence.jl")
include("bulk.jl")
include("holtslag.jl")
include("large.jl")
include("orad.jl")
include("grid.jl")
include("parameters.jl")

include("ForwardStep.jl")
include("SaveOutput.jl")

fname="out/test.nc"
ReadOutput(fname)

mtime=time

SST=TO[1]
ForwardStep(mtime)