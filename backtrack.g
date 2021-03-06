#############################################################################
##
##                             AutSemi package
##  ./backtrack.g
##                                                          Sergio Siccha
##
##  Copyright...
##
##  Compute the semigroup of partial isomorphisms via backtrack
##
#############################################################################
computePartialIsoSemigroup := function( args... )
  local graph, gens, allVertices, idempotents, G, countAdded,
  AddToG, root, accept,
  getChild, getNext, getParent, backtrack,
  timesSmallerGenSet,
  timesSizes,
  count;
  count := rec( accept := 0, reject := 0 );

  if Length( args ) = 1 then
    graph := args[1];
  elif Length( args ) = 2 then
    graph := args[1];
    gens := args[2];
  fi;

  timesSmallerGenSet := [];
  timesSizes := [];

  allVertices := Vertices( graph );
  # Initialise fullDistances of graph
  graph.fullDistances := computeDistances( graph );

  if IsBound( gens ) then
    G := InverseSemigroup( gens );
  else
    # If no generators were specified, G starts containing all idempotents
    idempotents := List( allVertices, x -> Difference( allVertices, [x] ) );
    idempotents := List( idempotents, x -> PartialPerm( x, x ) );
    G := InverseSemigroup( idempotents );
  fi;

  countAdded := 0;
  ##TODO better heuristic for when to compute smaller generating set
  AddToG := function( x )
    local gens, res, time;
    if not x in G then
      ## heuristic for when to compute smaller generating set
      if false then
      #if countAdded mod 5 = 4 then
        #gens := ShallowCopy( SmallInverseSemigroupGeneratingSet( G ) );
        Print( "Computing smaller GenSet...\c" );
        res := GET_REAL_TIME_OF_FUNCTION_CALL(
            SmallInverseSemigroupGeneratingSet,
            [G],
            rec( passResult := true )
        );
        time := QuoInt( res.time, 10^6 );
        Print( time, "s\n" );
        Add( timesSmallerGenSet, time );
        gens := ShallowCopy( res.result );
        G := InverseSemigroup( gens );
      fi;
      gens := ShallowCopy( GeneratorsOfInverseSemigroup( G ) );
      #G := InverseSemigroup( gens );
      ## ClosureInverseSemigroup reuses information about G
      G := ClosureInverseSemigroup( G, x );
      Add( gens, x );
      Print( "Computing size...\c" );
      time := GET_REAL_TIME_OF_FUNCTION_CALL( Size, [G] );
      time := QuoInt( time, 10^6 );
      Print( time, "s\n" );
      Add( timesSizes, time );
      Print("--- size: ", Size( G ), ", gens: ", Length(gens), " ---\n");
      countAdded := countAdded + 1;
    fi;
  end;

  root := function()
    return EmptyPartialPerm();
  end;

  accept := function( x )
    return IsPartialAut( graph, x );
  end;

  # Add smallest possible point to domain of x
  getChild := function( x )
    local maxDomain, newPoint, newImage;
    if allVertices = DomainOfPartialPerm( x ) then
      return fail;
    fi;
    # always only add bigger points to domain
    if x = EmptyPartialPerm() then
      maxDomain := 0;
    else
      maxDomain := Maximum( DomainOfPartialPerm( x ) );
    fi;
    if maxDomain = Maximum( allVertices ) then
      return fail;
    fi;
    newPoint := Minimum( Filtered( allVertices, x -> (x > maxDomain) ) );
    newImage := Minimum( Difference( allVertices, ImageListOfPartialPerm( x ) ) );
    x := JoinOfPartialPerms(
      x,
      PartialPerm( [ newPoint ], [ newImage ] )
    );
    return x;
  end;

  # Map the biggest point of domain to the next point
  getNext := function( x )
    local maxDomain, possibleImages, img, posImg, nextImg;
    # x already permutes all points, there is no 'next'
    if allVertices = DomainOfPartialPerm( x ) then
      return fail;
    fi;
    maxDomain := Maximum( DomainOfPartialPerm( x ) );
    img := maxDomain ^ x;
    # restrict to parent
    x := RestrictedPartialPerm(
      x,
      Difference( DomainOfPartialPerm( x ), [ maxDomain ] )
    );
    possibleImages := Difference( allVertices, ImageSetOfPartialPerm( x ) );
    posImg := PositionSet( possibleImages, img );
    # Image of maxDomain can be increased
    if posImg < Size( possibleImages ) then
      nextImg := possibleImages[ posImg+1 ];
    # Image of maxDomain can not be increased anymore
    else
      # Try to increase maxDomain
      if maxDomain < Maximum( allVertices ) then
        maxDomain := maxDomain + 1;
        nextImg := Minimum( possibleImages );
      else
        return fail;
      fi;
    fi;
    # Assign image of maxDomain to the next point
    x := JoinOfPartialPerms(
      x,
      PartialPerm( [ maxDomain ], [ nextImg ] )
    );
    return x;
  end;

  getParent := function( x )
    local maxDomain;
    maxDomain := Maximum( DomainOfPartialPerm( x ) );
    x := RestrictedPartialPerm(
      x,
      Difference( DomainOfPartialPerm( x ), [ maxDomain ] )
    );
    return x;
  end;

  backtrack := function( x )
    local child;
    # If x is not a partial automorphism, its parent can be added to the semigroup
    if not accept( x ) then
      AddToG( getParent( x ) );
      count.reject := count.reject+1;
      return;
    fi;
    count.accept := count.accept+1;
    # Add x, if it is a leaf, i.e. a permutation of the whole graph
    if DomainOfPartialPerm( x ) = allVertices then
      AddToG( x );
    fi;
    child := getChild( x );
    while not child = fail do
      backtrack( child );
      child := getNext( child );
    od;
  end;

  backtrack( root() );
  return rec(
    G := G,
    timesSizes := timesSizes,
    timesSmallerGenSet := timesSmallerGenSet,
    count := count
  );
end;
