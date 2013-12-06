MatLab
======
Graphs and computations are based on the MATLAB programme below
with psi(beta,x), taben(n,x), sybt1d(rad,lgsph,ngauss), akin(n,x),
9
and aleig(a,b,eps) [1] not listed here.
clc;
close all;
clear all;
x=-1:0.1:1;
[y1,y2]=PH(x);
figure(1)
plot(x,y1,’-’,x,y2,’--’);
legend(’exact’,’approximated’);
title(’P_H(\mu) Exact and Approximated’);
xlabel(’\mu’);
ylabel(’P_H(\mu)’);
psi=psi(0.3055,0)
S=[0.2,0,0.5,0,0.3];
V=[0.4,0.7,0.4,1.3,2.1];
ref=1000;
nGauss=6;
nuSf=[1.4,0,0,0,0];
Ss=[0,0,0,0,0];
slab=pij_slab(S,V,ref,’tran’)
cylA=pij_cyl(S,V,6,’void’)
cylB=pij_cyl(S,V,6,’white’)
ref=10^(-5);
keffA=keff(cylB,V,nuSf,Ss,ref)
Ss=[0.05,0,0.05,0,0.05];
keffB=keff(cylB,V,nuSf,Ss,ref)
sphA=pij_sph(S,V,6,’void’)
sphB=pij_sph(S,V,6,’white’)
The functions below were used.
function [exact, approximated]=PH(muv)
n=size(muv,2);
for i=1:n
mu=muv(i);
if ((mu>=-1)&(mu<0))
10
exact(i)=0;
approximated(i)=1/2+mu;
elseif ((mu>=0)&(mu<=1))
exact(i)=2*mu;
approximated(i)=1/2+mu;
end
end
function f=pij_slab(S,V,ref,bc)
n=length(V);
Lcell=sum(S.*V);
if strcmp(bc,’void’)
for i=1:n
f(i,i)=1/V(i)*Rii(i,0,1,S,V,Lcell);
for j=i+1:n
f(i,j)=1/V(i)*Rij(i,j,0,1,S,V,Lcell);
f(j,i)=V(i)/V(j)*f(i,j);
end
end
elseif strcmp(bc,’tran’)
for i=1:n
tmp=0;
for m=0:ref
tmp=tmp+Rii(i,m,1,S,V,Lcell)+Rii(i,m+1,-1,S,V,Lcell);
end
f(i,i)=1/V(i)*tmp;
for j=i+1:n
tmp=0;
for m=0:ref
tmp=tmp+Rij(i,j,m,1,S,V,Lcell)+Rij(i,j,m+1,-1,S,V,Lcell);
end
f(i,j)=1/V(i)*tmp;
f(j,i)=V(i)/V(j)*f(i,j);
end
end
else
error(’Unknown boundary condition!’);
end
function f=Rii(i,m,sign,S,V,Lcell)
if (S(i)~=0)
f=sign*V(i)*taben(2,m*Lcell)/S(i)-1/(S(i))^2 ...
*(taben(3,m*Lcell)-taben(3,m*Lcell+sign*L(i-1/2,i+1/2,S,V)));
elseif (m~=0)
f=taben(1,m*Lcell)*(V(i))^2/2;
11
else
f=Inf;
end
function f=L(x,y,S,V)
f=0;
for i=1:length(S)
if (i-1/2>=x & i+1/2<=y)
f=f+S(i)*V(i);
end
end
function val=Rij(i,j,m,sign,S,V,Lcell)
if (S(i)~=0 & S(j)~=0)
val=1/(2*S(i)*S(j))*(taben(3,m*Lcell+sign*L(i+1/2,j-1/2,S,V)) ...
-taben(3,m*Lcell+sign*L(i+1/2,j+1/2,S,V)) ...
+taben(3,m*Lcell+sign*L(i-1/2,j+1/2,S,V)) ...
-taben(3,m*Lcell+sign*L(i-1/2,j-1/2,S,V)));
elseif (S(i)==0 & S(j)~=0)
val=sign*V(i)/(2*S(j))*(taben(2,m*Lcell+sign*L(i-1/2,j-1/2,S,V)) ...
-taben(2,m*Lcell+sign*L(i-1/2,j+1/2,S,V)));
elseif (S(j)==0 & S(i)~=0)
val=sign*V(j)/(2*S(i))*(taben(2,m*Lcell+sign*L(i+1/2,j-1/2,S,V)) ...
-taben(2,m*Lcell+sign*L(i-1/2,j-1/2,S,V)));
else
val=V(i)*V(j)/2*taben(1,m*Lcell+sign*L(i-1/2,j-1/2,S,V));
end
function f=pij_cyl(S,V,nGauss,bc)
r(1)=sqrt(V(1)/pi);
for z=2:length(V)
r(z)=sqrt(V(z)/pi+r(z-1)^2);
end
track=sybt1d(r,false,nGauss);
if strcmp(bc,’void’)
f=Lij(’cyl’,track,S,V,nGauss);
elseif strcmp(bc,’white’)
tmp=Lij(’cyl’,track,S,V,nGauss);
PiS=1-tmp*S’;
A=2*pi*r(length(r));
pSi=4.*V’.*PiS/A;
PSS=1-pSi’*S’;
f=tmp+1/(1-PSS)*PiS*pSi’;
else
error(’Unknown boundary condition!’);
12
end
function f=Lij(geom,track,S,V,nGauss)
n=length(V);
for i=1:n
tmp1=0;
tmp2=0;
for k=1:nGauss
if i>1
for m=1:i-1
tmp1=tmp1+track{m,k}{1}*(2*Di(geom,i,m,k,track,S) ...
+Cij(geom,tau(i,i,m,k,track,S),i,i,m,k,track,S));
end
end
tmp2=tmp2+track{i,k}{1}*Di(geom,i,i,k,track,S);
end
f(i,i)=1/V(i)*(tmp1+tmp2);
for j=i+1:n
tmp1=0;
tmp2=0;
for k=1:nGauss
if i>1
for m=1:i-1
x=i-m+1;
tmp1=tmp1+track{m,k}{1} ...
*(Cij(geom,tau(i,j,m,k,track,S) ...
+tau(i,i,m,k,track,S)+S(i) ...
*track{m,k}{3}(x),i,j,m,k,track,S)+...
Cij(geom,tau(i,j,m,k,track,S),i,j,m,k,track,S));
end
end
tmp2=tmp2+track{i,k}{1} ...
*Cij(geom,tau(i,j,i,k,track,S),i,j,i,k,track,S);
end
f(i,j)=1/V(i)*(tmp1+tmp2);
f(j,i)=V(i)/V(j)*f(i,j);
end
end
function f=Di(geom,i,m,nGauss,track,S)
x=i-m+1;
if strcmp(geom,’cyl’)
if S(i)~=0
f=Lx(x,m,nGauss,track)/S(i)-1/S(i)^2 ...
*(akin(3,0)-akin(3,S(i)*Lx(x,m,nGauss,track)));
13
else
f=pi*Lx(x,m,nGauss,track)^2/4;
end
elseif strcmp(geom,’sph’)
if S(i)~=0
f=Lx(x,m,nGauss,track)/S(i)-1/S(i)^2 ...
*(1-exp(-(S(i)*Lx(x,m,nGauss,track))));
else
f=Lx(x,m,nGauss,track)^2/2;
end
else
error(’Undefined geometry!’);
end
function f=Lx(n,m,nGauss,track)
if n==1
f=2*track{m,nGauss}{3}(n);
else
f=track{m,nGauss}{3}(n);
end
function f=Cij(geom,tau,i,j,m,nGauss,track,S)
xi=i-m+1;
xj=j-m+1;
if strcmp(geom,’cyl’)
if i<=j && S(i)~=0 && S(j)~=0
f=1/(S(i)*S(j))*(akin(3,tau)-akin(3,tau+S(i) ...
*Lx(xi,m,nGauss,track)) ...
-akin(3,tau+S(j)*Lx(xj,m,nGauss,track)) ...
+akin(3,tau+S(i)*Lx(xi,m,nGauss,track) ...
+S(j)*Lx(xj,m,nGauss,track)));
elseif i<=j && S(i)==0 && S(j)~=0
f=Lx(xi,m,nGauss,track)/S(j) ...
*(akin(2,tau)-akin(2,tau+S(j)*Lx(xj,m,nGauss,track)));
elseif i<=j && S(j)==0 && S(i)~=0
f=Lx(xj,m,nGauss,track)/S(i) ...
*(akin(2,tau)-akin(2,tau+S(i)*Lx(xi,m,nGauss,track)));
elseif i<=j && S(j)==0 && S(i)==0
f=Lx(xi,m,nGauss,track)*Lx(xj,m,nGauss,track)*akin(1,tau);
end
elseif strcmp(geom,’sph’)
if i<=j && S(i)~=0 && S(j)~=0
f=1/(S(i)*S(j))*(exp(-tau)-exp(-(tau+S(i) ...
*Lx(xi,m,nGauss,track)))- ...
exp(-(tau+S(j)*Lx(xj,m,nGauss,track))) ...
14
+exp(-(tau+S(i)*Lx(xi,m,nGauss,track) ...
+S(j)*Lx(xj,m,nGauss,track))));
elseif i<=j && S(i)==0 && S(j)~=0
f=Lx(xi,m,nGauss,track)/S(j) ...
*(exp(-tau)-exp(-(tau+S(j)*Lx(xj,m,nGauss,track))));
elseif i<=j && S(j)==0 && S(i)~=0
f=Lx(xj,m,nGauss,track)/S(i) ...
*(exp(-tau)-exp(-(tau+S(i)*Lx(xi,m,nGauss,track))));
elseif i<=j && S(j)==0 && S(i)==0
f=Lx(xi,m,nGauss,track)*Lx(xj,m,nGauss,track)*exp(-tau);
end
else
error(’Undefined geometry!’);
end
function f=tau(i,j,m,nGauss,track,S)
f=0;
x=i-m+1;
if i==j
for n=1:x-1
f=f+S(m+n-1)*track{m,nGauss}{3}(n);
end
f=2*f;
elseif j>i
for n=1:j-i-1
f=f+S(i+n)*track{m,nGauss}{3}(x+n);
end
end
function f=keff(pij,V,nuSf,Ss,ref)
n=length(V);
for i=1:n
for j=1:n
SQ(i,j)=Ss(j)*V(j)*pij(j,i)/V(i);
FQ(i,j)=nuSf(j)*V(j)*pij(j,i)/V(i);
end
end
SQ=eye(n)-SQ;
tmp=aleig(SQ,FQ,eps);
f=1/tmp{2};
function f=pij_sph(S,V,nGauss,bc)
r(1)=(3*V(1)/(4*pi))^(1/3);
for z=2:length(V)
r(z)=(3*V(z)/(4*pi)+r(z-1)^3)^(1/3);
15
end
track=sybt1d(r,true,nGauss);
if strcmp(bc,’void’)
f=Lij(’sph’,track,S,V,nGauss);
elseif strcmp(bc,’white’)
tmp=Lij(’sph’,track,S,V,nGauss);
PiS=1-tmp*S’;
A=2*pi*r(length(r))^2;
pSi=4.*V’.*PiS/A;
PSS=1-pSi’*S’;
f=tmp+1/(1-PSS)*PiS*pSi’;
else
error(’Unknown boundary condition!’);
end
1
