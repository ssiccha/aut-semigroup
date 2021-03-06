## How to get reps for orbits

GetVertexNames := function( graph )
    return graph.names;
end;


partition := function( graph, G )
  local orb, stronglyConnectedComponents,
  repsVertices, repsPositions, reps,
  repsGraphs, repsGraphsConnected, repsConnected,
  sizes, sizesConnected;
  orb := Orb( G, Vertices( graph ), OnSets, rec( orbitgraph := true ) );
  stronglyConnectedComponents := OrbSCC( orb );
  repsPositions := List( stronglyConnectedComponents, x -> x[1] );
  reps := orb{ repsPositions };
  repsGraphs := List( reps, x -> InducedSubgraph( graph, x ) );
  repsGraphsConnected := Filtered( repsGraphs, IsConnectedGraph );
  repsConnected := List( repsGraphsConnected, GetVertexNames );
  sizes := Collected( List( reps, Length ) );
  sizesConnected := Collected( List( repsConnected, Length ) );

  return rec(
    orb := orb,
    scc := stronglyConnectedComponents,
    repsPositions := repsPositions,
    reps := reps,
    repsConnected := repsConnected,
    sizes := sizes,
    sizesConnected := sizesConnected
  );
end;

partition8 := function( subgraph )
  local G, orb;
  G := InverseSemigroup( autGens8 );
  orb := Orb( G, subgraph, OnSets, rec( orbitgraph := true ) );
  return rec( orb := orb, scc := OrbSCC( orb ) );
end;
