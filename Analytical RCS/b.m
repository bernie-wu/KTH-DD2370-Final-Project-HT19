function blm = b(l,ka)
RBj = @(n,x)sqrt(pi*x/2).*besselj(n+1/2,x);
RBh = @(n,x)sqrt(pi*x/2).*besselh(n+1/2,2,x);
bl = (-1i)^l*sqrt(pi*(2*l+1)/(l*(l+1)));
U = (ka.*RBj(l-1,ka))-(l*RBj(l,ka));
V = (ka.*RBh(l-1,ka))-(l*RBh(l,ka));
blm = -bl.*U./V;
end