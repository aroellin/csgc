
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
# A named vector with basic centred subgraph count statistics
  
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

  # Calculat esome intermediary matrices
  Yh.2p = Yh %*% Yh
  Yh.3p = Yh.2p %*% Yh
  Yh.4p = Yh.3p %*% Yh
  Yh.2 = Yh^2
  Yh.3 = Yh^3
  Yh.4 = Yh^4
  vr.2p = vr %*% vr
  vr.3p = vr.2p %*% vr
  vr.4p = vr.3p %*% vr
  vr.2 = vr^2
  vr.3 = vr^3
  vr.4 = vr^4
  
  #
  V.edge = sum(Yh[upper.tri(Yh)]) 
  vr.edge = sum(vr[upper.tri(vr)])
  #
  V.twostar = sum((Yh.2p)[upper.tri(Yh)]) 
  vr.twostar = sum((vr.2p)[upper.tri(vr)]) 
  #
  V.triangle = sum(diag(Yh.3p))/6 
  vr.triangle = sum(diag(vr.3p))/6 
  #
  V.4cycle = (trace(Yh.4p) - 4 * sum((Yh.2 %*% Yh.2)[upper.tri(Yh)]) +
              - 2*sum(Yh.4)/2 )/8
  vr.4cycle = (trace(vr.4p) +
              - 4 * sum((vr.2 %*% vr.2)[upper.tri(vr)]) - 2*sum(vr.4)/2 )/8
  #
  V.3path = sum((Yh.3p)[upper.tri(Yh)]) +
              + sum((Yh.3)[upper.tri(Yh)]) - sum((Yh.2 %*% Yh))
  vr.3path = sum((vr.3p)[upper.tri(vr)]) +
              + sum((vr.3)[upper.tri(vr)]) - sum((vr.2 %*% vr))
  #
  V.3star = (sum(rowSums(Yh)^3) - 3*sum(rowSums(Yh.2)*rowSums(Yh)) + 2*sum(rowSums(Yh.3)))/6
  vr.3star = (sum(rowSums(vr)^3) - 3*sum(rowSums(vr.2)*rowSums(vr)) + 2*sum(rowSums(vr.3)))/6

  z.edge = V.edge / sqrt(vr.edge)
  z.twostar = V.twostar / sqrt(vr.twostar)
  z.triangle = V.triangle / sqrt(vr.triangle)
  z.4cycle = V.4cycle / sqrt(vr.4cycle)
  z.3path = V.3path / sqrt(vr.3path)
  z.3star = V.3star / sqrt(vr.3star)
 
  V = c(V.edge, V.twostar, V.triangle, V.4cycle, V.3path, V.3star)
  names(V) = c("V.edge", "V.twostar", "V.triangle", "V.4cycle", "V.3path", "V.3star")
  vr = c(vr.edge, vr.twostar, vr.triangle, vr.4cycle, vr.3path, vr.3star)
  names(vr) = c("var.edge", "var.twostar", "var.triangle", "var.4cycle", "var.3path", "var.3star")
  z = V/sqrt(vr)
  names(z) = c("z.edge", "z.twostar", "z.triangle", "z.4cycle", "z.3path", "z.3star")
  ans = list(z=z, V=V, var=vr)
  return(ans)
}



