## Copyright (C) 2020 Fredy Vides
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## 
## [K1718,A,x17,x18,xx,yy,zz,f]=DiffCOVID19(17,18,1,eps,1);
##
## Example:
## [K1718,A,x0,x,xx,yy,zz,f]=DiffCOVID19(17,18,1,eps,1);

## Author: fredy <fredy@HPCLAB>
## Created: 2020-04-01

function [K,A,x0,x,xx,yy,zz,f]=DiffCOVID19(m,n,p,tol,graph)
m=m+1;
n=n+1;
pkg load io;
COVIDHist=xlsread ('COVID19History.xlsx');
%HNConnect=xlsread ('HNConnect0.xlsx');
HNConnect=xlsread ('HNConnect1.xlsx');
A=HNConnect (1:18,1:18);
[M,N]=size(A);
E=eye(M,N);
if nargin<=4
	graph=1;
end
x0=COVIDHist (1:18,m);
x=COVIDHist (1:18,n);
K=diag(x0>0)*A*diag(1./sum(A));
K=E-K;
d=K*x0;
d=-(x-x0)'*d/(d'*d);
K=E-d*K;
[xx,yy]=meshgrid(-1.5:3/30:1.5);
f=[];
r=.5;
	z1=(r*exp(2*pi*i*(0:6)/7)).';
	z2=(2.0*r*exp(20*pi*i*(0:8)/(9*21))).';
	z3=2.4*r*exp((pi+.1)*i/4);
	xy=zeros(M,2);
	xy([15 18 4 12 17 2 7],:)=[real(z1),imag(z1)];
	xy([9 3 1 6 16 5 14 13 10],:)=[real(z2),imag(z2)];
	xy(11,:)=[real(z3),imag(z3)];
for k=1:18
	z=sqrt((xx-xy(k,1)).^2+(yy-xy(k,2)).^2);
	m=min(min(z));
	f=[f;find(z==m)];
end
  [X,Y]=meshgrid (-1.5:3/150:1.5);
  zz=zeros(size(xx));
  zz(f)=Mxvec(K,x0,p);
  Z=interp2(xx,yy,zz,X,Y,'spline');
if graph==1
  figure;
  contour(X,Y,Z,64);
  colormap summer;
  hold on;
  gplot (A,[xx(f) yy(f)],'k-');
  plot(xx(f),yy(f),'k.','markersize',20,xx(f(8)),...
  yy(f(8)),'r.','markersize',20,xx(f(6)),yy(f(6)),...
  'b.','markersize',20);
  hold off;
  axis off;
  axis square;
end
end
function y=Mxvec(A,b,n)
	y=b;	
	for k=1:n
		y=A*y;
	end
end