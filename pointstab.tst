
BruteForcePointStabilizer := function( S, pnt )
    local l;
    l := AsList( S );
    return SortedList( Filtered( l, x -> pnt^x = pnt ) );
end;

test := function( S, pnt )
    local res, bf;
    if not pnt in MovedPoints(S) then
        return true;
    fi;
    res := PointStabilizer( S, pnt );
    bf := BruteForcePointStabilizer( S, pnt );
    res := AsSortedList( res );
    Print(Size(S), ",");
    return bf = res;
end;
