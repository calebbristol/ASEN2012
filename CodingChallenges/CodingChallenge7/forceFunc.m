function F = forceFunc(t,V,Cd,rho,A,mass,g)
D = Cd * 0.5 * rho * V^2 * A;
F = D - mass*g;
end

