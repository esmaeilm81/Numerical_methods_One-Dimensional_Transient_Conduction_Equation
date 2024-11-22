clc;
clear;
close all;

h=0.015625;
N=(1/h)+1;
alpha=1;
v=10;
t=v*h^2/alpha;
tol=1e-9;
maxit=N;

B=zeros(N,N);
phi_p=zeros(N,1);
phi_p(1)=1;
phi_p(N)=0;
for j=1:N
   x(j)=(j-1)*h;
end
B(1,1)=1; B(N,N)=1;
k=2;
for j=2:N-1
   B(k,j-1)=-v;
   B(k,j)=1+2*v;
   B(k,j+1)=-v;
   k=k+1;
end
B=sparse(B);
for i=1:10000000
    phi=gmres(B,phi_p,[],tol,maxit);
    error(i)=max(abs(phi-phi_p));
    if error(i)<2e-09
        break;
    end
    phi_p=phi;
end
figure();
plot(x',phi);
for i=1:length(error)
    T(i)=i*t;
end
figure();
plot(T,error);
