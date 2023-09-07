# This is a companion file of the article "Computations of Gromov--Witten invariants of toric varieties" by Giosuè Muratore.
# In this file are reproduced all computations in that article.
# In order to execute this file, you need the Julia package "ToricAtiyahBott".

# Section 4.1
P3 = projective_space(NormalToricVariety, 3); # the space P^3.
h = cohomology_class(toric_line_bundle(P3, [1]));
P = ev(1, h^3)*ev(2, h^3)*ev(3, h^3)*ev(4, h^2)*ev(5, h^2);
IntegrateAB(P3, 2*h^2, 5, P);

X = domain(blow_up(P3, ideal(gens(cox_ring(P3))[1:3])));
mg = moment_graph(X);
H = mg[4,5];
E = mg[1,2];
d = 1; # this can be any non negative integer
e = -1; # this can be any non positive integer
beta = d*H+e*E;
P = prod(i -> ev(i, a_point(X)), 1:(2*d+e)); # we require 2d+e>0
IntegrateAB(X, beta, 2*d+e, P);

l = 4*h;
P = 1//4*(ev(1,l)^3+3*ev(1,l)^2*Psi(1)+2*ev(1,l)*Psi(2))*ev(1,h)^2;
IntegrateAB(P3, h^2, 1, P);

# Section 4.2
e1 = [1,0,0]; e2 = [0,1,0]; e3 = [0,0,1];
ray_gens = [e1, -e1, e2, e3, -e2-e3-2*e1];
max_cones = [[1,3,4],[1,3,5],[1,4,5],[2,3,4],[2,3,5],[2,4,5]];
X=normal_toric_variety(ray_gens,max_cones;non_redundant=true);
Z = [cohomology_class(X, Z) for Z in gens(cohomology_ring(X))];
mg = moment_graph(X);
lambda1 = mg[1,4];
lambda2 = mg[4,5];
P = ev(1, Z[3])*ev(2, Z[3])*ev(3, Z[4]*Z[5]);
IntegrateAB(X, lambda2, 3, P);
P = ev(1, Z[1])*ev(2, Z[1])*ev(3, Z[2]*Z[3]*Z[4]);
IntegrateAB(X, lambda1, 3, P);

# Section 4.3
X = proj(toric_line_bundle(P3, [0]),toric_line_bundle(P3, [1]));
mg = moment_graph(X);
lambda1 = mg[1,2];
lambda2 = mg[4,8];
integrate(cohomology_class(anticanonical_divisor(X))*lambda1)
integrate(cohomology_class(anticanonical_divisor(X))*lambda2)
Z = [cohomology_class(X, Z) for Z in gens(cohomology_ring(X))];
P = ev(1, Z[4])*ev(2, Z[5])*ev(3, a_point(X));
IntegrateAB(X, lambda1, 3, P);

# Section 4.4
M = anticanonical_bundle(X);
P = push_ev(M);
IntegrateAB(X, lambda1, 1, Psi(1)*P);
IntegrateAB(X, lambda2, 0, P);
