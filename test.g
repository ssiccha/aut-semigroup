Read("read.g");
Read("examples.g");
tmp := constructMeshAndAutSubsemi( 4 );
G := computePartialIsoSemigroup( tmp.mesh, tmp.autGens );;
