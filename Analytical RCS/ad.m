function alm = ad(l,ka)
RBj = @(n,x)sqrt(pi*x/2).*besselj(n+1/2,x);
RBh = @(n,x)sqrt(pi*x/2).*besselh(n+1/2,2,x);
global eta lambda;
k1 = sqrt(lambda);
kb = ka/eta;
k1a = k1*ka;
k1b = k1a/eta;
al = (-1i)^l*sqrt(pi*(2*l+1)/(l*(l+1)));
A1 = RBj(l,k1a); A2 = RBj(l,k1b); A3 = ((k1b.*RBj(l-1,k1b))-(l*RBj(l,k1b)))/k1;
B1 = RBh(l,k1a); B2 = RBh(l,k1b); B3 = ((k1b.*RBh(l-1,k1b))-(l*RBh(l,k1b)))/k1; 
C2 = RBj(l,kb); C3 = (kb.*RBj(l-1,kb))-(l*RBj(l,kb));
D2 = RBh(l,kb); D3 = (kb.*RBh(l-1,kb))-(l*RBh(l,kb));
num = (A3.*B1.*C2)-(A1.*B3.*C2)-(A2.*B1.*C3)+(A1.*B2.*C3);
den = (A3.*B1.*D2)-(A1.*B3.*D2)-(A2.*B1.*D3)+(A1.*B2.*D3);
alm = -al*num./den;
end