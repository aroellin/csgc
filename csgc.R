
# Return trace of a matrix
trace = function(A) return(sum(diag(A)))

# Calculate basic centred subgraph counts
csgc = function(A, P) {
  #
  # INPUT
  # A     = adjacency matrix
  # P     = hypothesised connection probabilities
  # A and P need to have the same dimensions
  #
  # OUTPUT
  # A list with named vectors:
  # z:   normalised centred subgraph counts
  # V:   unnormalised centred subgraph counts
  # var: variances of V
  #
  # NOTE
  # csgc(A,0*A) will return the raw subgraph counts in V
  
  # Clean up A and P
  diag(A) = 0 # remove loops
  diag(P) = 0
  # ignore those entries where P>1
  w = which(P>1,arr.ind=TRUE)
  P[w] = 1
  
  # Calculate variance matrix
  vr = P*(1-P)
  # Calculate centred indicators
  Yh = A - P
  # Remove those where P>1
  Yh[w] = 0
  
  # Calculate some intermediary matrices
  Yh.2p = Yh %*% Yh
  Yh.3p = Yh.2p %*% Yh
  Yh.4p = Yh.3p %*% Yh
  Yh.5p = Yh.4p %*% Yh
  Yh.2 = Yh^2
  Yh.3 = Yh^3
  Yh.4 = Yh^4
  Yh.5 = Yh^5
  vr.2p = vr %*% vr
  vr.3p = vr.2p %*% vr
  vr.4p = vr.3p %*% vr
  vr.5p = vr.4p %*% vr
  vr.2 = vr^2
  vr.3 = vr^3
  vr.4 = vr^4
  vr.5 = vr^5
  
  #
  V.edge = sum(Yh[upper.tri(Yh)]) 
  vr.edge = sum(vr[upper.tri(vr)])
  #
  V.twostar = sum((Yh.2p)[upper.tri(Yh)]) 
  vr.twostar = sum((vr.2p)[upper.tri(vr)]) 
  #
  V.triangle = trace(Yh.3p)/6 
  vr.triangle = trace(vr.3p)/6 
  #
  V.4cycle = (trace(Yh.4p) - 4 * sum((Yh.2 %*% Yh.2)[upper.tri(Yh)]) - sum(Yh.4) )/8
  vr.4cycle = (trace(vr.4p) - 4 * sum((vr.2 %*% vr.2)[upper.tri(vr)]) - sum(vr.4) )/8
  #
  V.3path = sum((Yh.3p)[upper.tri(Yh)]) + sum((Yh.3)[upper.tri(Yh)]) - sum((Yh.2 %*% Yh))
  vr.3path = sum((vr.3p)[upper.tri(vr)]) + sum((vr.3)[upper.tri(vr)]) - sum((vr.2 %*% vr))
  #
  V.3star = (sum(rowSums(Yh)^3) - 3*sum(Yh.2 %*% Yh) + 2*sum(Yh.3))/6
  vr.3star = (sum(rowSums(vr)^3) - 3*sum(vr.2 %*% vr) + 2*sum((vr.3)))/6
  #
  V.triangleappendix = (sum(diag(Yh.3p) %*% Yh) - 2*trace(Yh.2p %*% Yh.2))/2
  vr.triangleappendix = (sum(diag(vr.3p) %*% vr) - 2*trace(vr.2p %*% vr.2))/2
  #
  V.twotriangle = (sum(Yh.2p * Yh.2p * Yh) - trace(Yh.2 %*% Yh.2 %*% Yh))/4
  vr.twotriangle = (sum(vr.2p * vr.2p * vr) - trace(vr.2 %*% vr.2 %*% vr))/4
  #
  V.fivecycle = (trace(Yh.5p) - 5*sum(diag(Yh.3p)%*%Yh.2) + 5*trace(Yh.3%*%Yh.2p))/10
  vr.fivecycle = (trace(vr.5p) - 5*sum(diag(vr.3p)%*%vr.2) + 5*trace(vr.3%*%vr.2p))/10
  #
  V.fourpath = (2*sum((Yh.4p)[upper.tri(Yh)]) - 2*sum(Yh.2 %*% Yh.2p) +
                  2*sum(Yh.3 %*% Yh) + 3*sum(Yh.2 %*% Yh.2) + 3*trace(Yh.2 %*% Yh.2p) -
                  2*sum(Yh.4) - sum(rowSums(Yh)^2 * rowSums(Yh.2)) - 2*sum(diag(Yh.3p) %*% Yh))/2
  vr.fourpath = (2*sum((vr.4p)[upper.tri(vr)]) - 2*sum(vr.2 %*% vr.2p) +
                   2*sum(vr.3 %*% vr) + 3*sum(vr.2 %*% vr.2) + 3*trace(vr.2 %*% vr.2p) -
                   2*sum(vr.4) - sum(rowSums(vr)^2 * rowSums(vr.2)) - 2*sum(diag(vr.3p) %*% vr))/2
  
  z.edge = V.edge / sqrt(vr.edge)
  z.twostar = V.twostar / sqrt(vr.twostar)
  z.triangle = V.triangle / sqrt(vr.triangle)
  z.4cycle = V.4cycle / sqrt(vr.4cycle)
  z.3path = V.3path / sqrt(vr.3path)
  z.3star = V.3star / sqrt(vr.3star)
  z.triangleappendix = V.triangleappendix / sqrt(vr.triangleappendix)
  z.twotriangle = V.twotriangle / sqrt(vr.twotriangle)
  z.fivecycle = V.fivecycle / sqrt(vr.fivecycle)
  z.fourpath = V.fourpath / sqrt(vr.fourpath)
  
  V = c(V.edge, V.twostar, V.triangle, V.4cycle, V.3path, V.3star, 
        V.triangleappendix, V.twotriangle, V.fivecycle, V.fourpath)
  names(V) = c("V.edge", "V.twostar", "V.triangle", "V.4cycle", "V.3path", "V.3star",
               "V.triangleappendix", "V.twotriangle", "V.fivecycle", "V.fourpath")
  vr = c(vr.edge, vr.twostar, vr.triangle, vr.4cycle, vr.3path, vr.3star,
         vr.triangleappendix, vr.twotriangle, vr.fivecycle, vr.fourpath)
  names(vr) = c("var.edge", "var.twostar", "var.triangle", "var.4cycle", "var.3path", "var.3star",
                "var.triangleappendix", "var.twotriangle", "var.fivecycle", "vr.fourpath")
  z = V/sqrt(vr)
  names(z) = c("z.edge", "z.twostar", "z.triangle", "z.4cycle", "z.3path", "z.3star",
               "z.triangleappendix", "z.twotriangle", "z.fivecycle", "z.fourpath")
  ans = list(z=z, V=V, var=vr)
  return(ans)
}



