Read("read.g");
res := constructMeshAndAutSubsemi( 4 );;
G := computePartialIsoSemigroup( res.mesh, res.autGens );; time;
