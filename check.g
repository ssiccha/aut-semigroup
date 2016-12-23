###############################
# function IsPartialAut
# Input:
#   x - PartialPerm
#
# Output:
#   whether x is a partial automorphism of graph
###############################
IsPartialAut := function( graph, x )
  if x = EmptyPartialPerm() then
    return true;
  fi;
  return
    lookupDistances( graph, DomainOfPartialPerm(x) )
    =
    lookupDistances( graph, ImageListOfPartialPerm(x) );
end;
