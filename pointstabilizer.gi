PointStabilizer := function( S, POINT )
    local orb, sccLookup, setsWithPoint, sccsOfSetsWithPoint, toSCCRep, fromSCCRep, stabGens, c, d, stab, newStabGens, dom, i, pointStabilizer;
    orb := LambdaOrb( S );
    sccLookup := OrbSCCLookup( orb );
    setsWithPoint := Filtered( [1..Size(orb)], i -> POINT in orb[i] );
    sccLookup := sccLookup{ setsWithPoint };
    sccReps := List( sccLookup, i -> orb[ OrbSCC(orb)[i][1] ] );
    tmp := List(
        [1..Size(setsWithPoint)],
        i -> LambdaOrbMult(LambdaOrb(S), sccLookup[i], setsWithPoint[i] )
    );
    fromSCCRep := List( tmp, x -> x[1] );
    toSCCRep := List( tmp, x -> x[2] );
    stabGens := [];
    for i in [ 1 .. Length(setsWithPoint) ] do
        to := toSCCRep[i];
        from := fromSCCRep[i];
        # Let x be the set containing POINT. To compute Stab_POINT in the Schutzenberger
        # group of x, we need to work with the set representing the SCC of x.
        stab := Stabilizer( LambdaOrbSchutzGp(orb, sccLookup[i]), POINT^to );
        # Convert the permutations in stab to PartialPerms
        newStabGens := GeneratorsOfGroup( stab );
        dom := orb[OrbSCC(orb)[sccLookup[i]][1]];
        newStabGens := List( newStabGens, x -> AsPartialPerm( x, dom ) );
        # Conjugate them back into the Schutzenberger group of x
        newStabGens := List( newStabGens, x -> x^from );
        Add( stabGens, newStabGens );
    od;
    stabGens := Concatenation( stabGens );
    if IsInverseSemigroup(S) then
        pointStabilizer := InverseSemigroup( stabGens[1] );
        pointStabilizer := ClosureInverseSemigroup( pointStabilizer, stabGens );
    else
        pointStabilizer := Semigroup( stabGens[1] );
        pointStabilizer := ClosureSemigroup( pointStabilizer, stabGens );
    fi;
    return pointStabilizer;
end;
