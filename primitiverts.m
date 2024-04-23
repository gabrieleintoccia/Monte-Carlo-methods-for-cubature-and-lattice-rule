function g = primitiverts(p)
a=primes(p);
prm=a(length(a));
if p~=prm
    error ('n non è numero primo')
end
primef=unique(factor(p-1));
g = 2; 
i =1;
while i <= length(primef)
    if 	powermod(g, (p-1)/primef(i), p) == 1
    g = g + 1; i =0; 
    end
    i=i+1;
end