function [z, e2] = fastrank1(n, s_max, omega, gamma, beta)
a=primes(n);
prm=a(length(a));
if n~=prm
    error ('n non è numero primo')
end
z=zeros(s_max,1);
e2=zeros(s_max,1);
m = (n-1)/2;
E2=zeros(m,1);
cumbeta=cumprod(beta);

g = primitiverts(n);
perm=zeros(m,1);
perm(1)=1;
for j=1:m-1
    perm(j+1)=mod(perm(j)*g,n);
end
perm=min(n-perm,perm);
psi=omega(perm/n);
psi0=omega(0);
fft_psi=fft(psi);

q=ones(m,1);
q0=1;
for s=1:s_max
    E2=ifft(fft_psi.*fft(q));
    E2=real(E2);
    [min_E2,w(s)]=min(E2);
    z(s)=perm(w(s));
    e2(s)=-cumbeta(s)+(beta(s)*(q0+2*sum(q))+gamma(s)*(psi0*q0+2*min_E2))/n;
    q = (beta(s) + gamma(s) * psi([w(s):-1:1 m:-1:w(s)+1])) .* q;
    q0 = (beta(s) + gamma(s) * psi0)*q0;
    fprintf('s=%4d, z=%6d, e2=%.4e, e=%.4e\n',s,z(s),e2(s), sqrt(e2(s)));
end
