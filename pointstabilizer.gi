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
    for i in [ 1 .. Length(setsWithPoint) ] do
        for j in [ 1 .. Length(setsWithPoint) ] do
            if i neq j and sccLookup[i] eq sccLookup[j] then
                # For two different sets in the same SCC take semigroup elements
                # mapping one to the other via the SCC's representative.
                # Then check whether POINT^to and PreImage( from, POINT ) can be
                # mapped to each other in the representatives Schutzenberger
                # group.
                to := toSCCRep[i];
                from := fromSCCRep[j];
                schutzGp := LambdaOrbSchutzGp( orb, sccLookup[i] );
                if PreImagePartialPerm( from, POINT )
                in Orbit( schutzGp, POINT^to ) then
                    # Handle schutzGp not being transitive by switching
                    # to its constituent action on the orbit of POINT^to.
                    phi := ActionHomomorphism(
                        schutzGp,
                        Orbit( schutzGp, POINT^to ) );
                    stabChain := StabChainMutable( Image( phi ) );
                    # ActionHomomorphism renames the points via phi!.conperm
                    a := (POINT^to)^phi!.conperm;
                    b := (PreImagePartialPerm( from, POINT ))^phi!.conperm;
                    pi1 := InverseRepresentative( stabChain, a );
                    pi2 := InverseRepresentative( stabChain, b );
                    pi := pi1 / pi2;
                    pi := PreImagesRepresentative( phi, pi );
                    # Make this a partial perm again.
                    dom := orb[OrbSCC(orb)[sccLookup[i]][1]];
                    pi := AsPartialPerm( pi, dom );
                    Add( stabGens, to * pi * from );
                fi;
            fi;
        od;
    od;
    if IsInverseSemigroup(S) then
        pointStabilizer := InverseSemigroup( stabGens[1] );
        pointStabilizer := ClosureInverseSemigroup( pointStabilizer, stabGens );
    else
        pointStabilizer := Semigroup( stabGens[1] );
        pointStabilizer := ClosureSemigroup( pointStabilizer, stabGens );
    fi;
    return pointStabilizer;
end;
