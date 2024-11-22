clc;
clear;
close all;

h=0.015625;
N=(1/h)+1;
alpha=1;
v=0.25;
t=v*h^2/alpha;

A=zeros(N,N);
phi_p=zeros(N,1);
phi_p(1)=1;
phi_p(N)=0;
for j=1:N
   x(j)=(j-1)*h;
end
A(1,1)=1; A(N,N)=1;
k=2;
for j=2:N-1
   A(k,j-1)=v;
   A(k,j)=1-2*v;
   A(k,j+1)=v;
   k=k+1;
end
A=sparse(A);
for i=1:10000000
    phi=A*phi_p+B;
    error(i)=max(abs(phi-phi_p));
    if error(i)<2e-09
        break;
    elseif error(i)>100000
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


