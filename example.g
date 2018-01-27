Read("read.g");
Read("./construct-examples.g");
res := constructMeshAndAutSubsemi( 4 );;
G := computeAutSemigroup( res.mesh, res.autGens );; time;
