e1 = [1,0,0]; e2 = [0,1,0]; e3 = [0,0,1];
ray_gens = [e1, -e1, e2, e3, -e2-e3-2*e1];
max_cones = [[1,3,4],[1,3,5],[1,4,5],[2,3,4],[2,3,5],[2,4,5]];
X=normal_toric_variety(ray_gens,max_cones;non_redundant=true);

Z = [cohomology_class(X, Z) for Z in gens(cohomology_ring(X))];
mg = moment_graph(X);
lambda1 = mg[1,4];
lambda2 = mg[4,5];


#### lambda2 invariants
P = ev(1, Z[3])*ev(2, Z[3])*ev(3, Z[4]*Z[5]);
IntegrateAB(X, lambda2, 3, P); #-1


#### lambda1 invariants
P = ev(1, Z[1])*ev(2, Z[1])*ev(3, Z[2]*Z[3]*Z[4]);
Q = ev(1, Z[1])*ev(2, Z[1]*Z[4])*ev(3, Z[1]*Z[3]);
IntegrateAB(X, lambda1, 3, [P,Q]); #1;1


#### 2*lambda2 invariants
P = ev(1, Z[3])*ev(2, Z[3]*Z[4])*ev(3, Z[4]*Z[5]);
IntegrateAB(X, 2*lambda2, 3, P); #-2


#### 3*lambda2 invariants
P = ev(1, Z[3]*Z[4])*ev(2, Z[4]*Z[5])*ev(3, Z[3]*Z[5]);
IntegrateAB(X, 3*lambda2, 3, P); #-4


#### lambda1+lambda2 invariants
P = ev(1, Z[1])*ev(2, Z[4]*Z[5])*ev(3, Z[1]*Z[3]*Z[4]);
IntegrateAB(X, lambda1+lambda2, 3, P); #1


#### lambda1+2*lambda2 invariants
P = ev(1, Z[4]*Z[5])*ev(2, Z[3]*Z[5])*ev(3, Z[1]*Z[3]*Z[4]);
IntegrateAB(X, lambda1+2*lambda2, 3, P); #1


#### 2*lambda1+lambda2 invariants
P = ev(1, Z[2]*Z[2])*ev(2, Z[1]*Z[3]*Z[4])*ev(3, Z[1]*Z[4]*Z[5]);
Q = ev(1, Z[1]*Z[3])*ev(2, Z[1]*Z[3]*Z[4])*ev(3, Z[1]*Z[4]*Z[5]);
IntegrateAB(X, 2*lambda1+lambda2, 3, [P,Q]); #-2; 1