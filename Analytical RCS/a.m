function alm = a(l,ka)
RBj = @(n,x)sqrt(pi*x/2).*besselj(n+1/2,x);
RBh = @(n,x)sqrt(pi*x/2).*besselh(n+1/2,2,x);
al = (-1i)^l*sqrt(pi*(2*l+1)/(l*(l+1)));
A = RBj(l,ka);
B = RBh(l,ka);
alm = -al.*A./B;
end