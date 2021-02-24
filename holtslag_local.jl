function holtslag_local(Ri)
  if Ri>0
    y=1/(1 + 10 * Ri * (1 + 8 * Ri))
  else
    y=(1 -18 * Ri)^(0.5)
  end
  return y
end