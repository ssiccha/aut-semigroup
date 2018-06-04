Read("read.g");
Read("./construct-examples.g");
res := constructMeshAndAutSubsemi( 4 );;
G := computePartialIsoSemigroup( res.mesh, res.autGens );; time;
