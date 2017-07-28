Read("read.g");
Read("examples.g");
tmp := constructMeshAndAutSubsemi( 4 );
G := computeAutSemigroup( tmp.mesh, tmp.autGens );;
