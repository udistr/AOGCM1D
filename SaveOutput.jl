using DataStructures
using NCDatasets
using Dates

function SaveOutput(fname)
  Dataset(fname,"c",attrib = OrderedDict("title" => "test")) do ds
    defDim(ds,"Alon",1)
    defDim(ds,"Alat",1)
    defDim(ds,"Alev",vleva)
    defDim(ds,"Olon",1)
    defDim(ds,"Olat",1)
    defDim(ds,"Olev",vlevo)
    defDim(ds,"time",1)
    defVar(ds,"time",Dates.datetime2epochms(mtime)/1000/3600/24,("time",), 
      attrib = OrderedDict("units" => "days since 0000-01-01 00:00:00" ))
    defVar(ds,"Alon",lon,("Alon",), 
      attrib = OrderedDict("units" => "degree"))
    defVar(ds,"Alat",lon,("Alat",), 
      attrib = OrderedDict("units" => "degree"))
    defVar(ds,"Alev",ZAC,("Alev",), 
      attrib = OrderedDict("units" => "m"))
    defVar(ds,"Olon",lon,("Olon",), 
      attrib = OrderedDict("units" => "degree"))
    defVar(ds,"Olat",lon,("Olat",), 
      attrib = OrderedDict("units" => "degree"))
    defVar(ds,"Olev",ZOC,("Olev",), 
      attrib = OrderedDict("units" => "m"))
    defVar(ds,"THETA",ΘA,("Alev",), 
      attrib = OrderedDict("units" => "degree Kelvin"))
    defVar(ds,"TA",TA,("Alev",), 
      attrib = OrderedDict("units" => "degree Celsius"))
    defVar(ds,"qA",qA,("Alev",), 
      attrib = OrderedDict("units" => "kg/kg"))
    defVar(ds,"UA",UA,("Alev",), 
      attrib = OrderedDict("units" => "m/s"))
    defVar(ds,"rhoA",ρA,("Alev",), 
      attrib = OrderedDict("units" => "kg/m^3"))
    defVar(ds,"TO",TO,("Olev",), 
      attrib = OrderedDict("units" => "degree Celsius"))
    defVar(ds,"UO",UO,("Olev",), 
      attrib = OrderedDict("units" => "m/s"))
    defVar(ds,"rhoO",ρO,("Olev",), 
      attrib = OrderedDict("units" => "kg/kg"))
  end
end

function ReadOutput(fname)
  ds = NCDataset(fname)
  global lat=ds["Alat"][1]
  global lon=ds["Alon"][1]
  global time=ds["time"][1]
  global ZAC=ds["Alev"][:]
  global ZOC=ds["Olev"][:]
  global ΘA=ds["THETA"][:]
  global TA=ds["TA"][:]
  global qA=ds["qA"][:]
  global UA=ds["UA"][:]
  global ρA=ds["rhoA"][:]
  global TO=ds["TO"][:]
  global UO=ds["UO"][:]
  global ρO=ds["rhoO"][:]
end