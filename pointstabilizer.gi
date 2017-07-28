PointStabilizer := function( S, POINT )
    local setsOrb, sccLookup, setsWithPoint, indicesOfSCCs, sccsOfSetsWithPoint, canonizers, lambdaInversesOfCanonizers, stabGens, c, d, stab, newStabGens, dom, i, pointStabilizer;
    setsOrb := LambdaOrb( S );
    sccLookup := OrbSCCLookup( setsOrb );
    setsWithPoint := Filtered( [1..Size(setsOrb)], i -> POINT in setsOrb[i] );
    indicesOfSCCs := sccLookup{ setsWithPoint };
    sccsOfSetsWithPoint := sccLookup{ setsWithPoint };
    canonizers := List(
        [1..Size(setsWithPoint)],
        i -> LambdaOrbMult(LambdaOrb(S), indicesOfSCCs[i], setsWithPoint[i] )
    );
    lambdaInversesOfCanonizers := List( canonizers, x -> x[1] );
    canonizers := List( canonizers, x -> x[2] );
    stabGens := [];
    for i in [ 1 .. Length(setsWithPoint) ] do
        c := canonizers[i];
        d := lambdaInversesOfCanonizers[i];
        # Let x be the set containing POINT. To compute Stab_POINT in the Schutzenberger
        # group of x, we need to work with the set representing the SCC of x.
        stab := Stabilizer( LambdaOrbSchutzGp(setsOrb, indicesOfSCCs[i]), POINT^c );
        # Convert the permutations in stab to PartialPerms
        newStabGens := GeneratorsOfGroup( stab );
        dom := setsOrb[OrbSCC(setsOrb)[indicesOfSCCs[i]][1]];
        newStabGens := List( newStabGens, x -> PartialPerm( dom, OnTuples( dom, x ) ) );
        # Conjugate them back into the Schutzenberger group of x
        newStabGens := List( newStabGens, x -> x^d );
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
