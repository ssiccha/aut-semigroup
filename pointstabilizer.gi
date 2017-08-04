PointStabilizer := function( S, POINT )
    local orb, sccLookup, setsWithPoint, sccReps, tmp, fromSCCRep, toSCCRep, stabGens, to, from, stab, newStabGens, dom, schutzGp, phi, stabChain, a, b, pi1, pi2, pi, pointStabilizer, i, j;
    if ForAll( GeneratorsOfInverseSemigroup(S), x -> POINT^x = 0 ) then
        Error( "S must be defined on POINT!" );
    fi;
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
        # Let x be the set containing POINT. To compute Stab_POINT in the
        # Schutzenberger group of x, we need to work with the set representing
        # the SCC of x.
        schutzGp := LambdaOrbSchutzGp(orb, sccLookup[i]);
        stab := Stabilizer( schutzGp, POINT^to );
        # Convert the permutations in stab to PartialPerms
        newStabGens := GeneratorsOfGroup( stab );
        dom := orb[OrbSCC(orb)[sccLookup[i]][1]];
        newStabGens := List( newStabGens, x -> AsPartialPerm( x, dom ) );
        # If stab is trivial newStabGens will be empty. Add the identity pperm
        if IsEmpty( newStabGens ) then
            newStabGens := [ AsPartialPerm( (), dom ) ];
        fi;
        # Conjugate them back into the Schutzenberger group of x
        newStabGens := List( newStabGens, x -> x^(to^-1) );
        Add( stabGens, newStabGens );
    od;
    stabGens := Concatenation( stabGens );
    for i in [ 1 .. Length(setsWithPoint) ] do
        for j in [ 1 .. Length(setsWithPoint) ] do
            if not i = j and sccLookup[i] = sccLookup[j] then
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
                    # Handle schutzGp acting trivially on POINT^to
                    if IsTrivial( Image( phi ) ) then
                        Add( stabGens, to * from );
                        continue;
                    fi;
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
    # Try to keep the generating set small by using ClosureInverseSemigroup
    pointStabilizer := InverseSemigroup( stabGens[1] );
    for i in [ 1 .. QuoInt( Length( stabGens ), 3 ) ] do
        pointStabilizer := ClosureInverseSemigroup(
            pointStabilizer,
            stabGens{[3*i-2 .. 3*i]}
        );
    od;
    modulus := Length( stabGens ) mod 3;
    if not modulus = 0 then
        pointStabilizer := ClosureInverseSemigroup(
            pointStabilizer,
            stabGens{[ Length( stabGens ) - modulus + 1, Length( stabGens ) ]}
        );
    fi;
    return pointStabilizer;
end;
