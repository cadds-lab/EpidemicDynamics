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
## function [K,T,x0,x]=UACCorrector(m,n,tol,graph)
##
## Example:
## [K,T01,x0,x1]=UACCorrector(0,1,eps);
## [K,T12,x1,x2]=UACCorrector(1,2,eps);
## [K,T23,x2,x3]=UACCorrector(2,3,eps);
## [K,T03,x0,x3]=UACCorrector(0,3,eps);
## norm(x3-T03*x0,1)+norm((T23*T12*T01-T03)*x0,1)

## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-17

function [K,T,x0,x]=UACCorrector(m,n,tol,graph)
m=m+1;
n=n+1;
pkg load io;
COVIDHist=xlsread ('COVID19HistoryJoint.xlsx');
%HNConnect=xlsread ('HNConnect0.xlsx');
HNConnect=xlsread ('HNConnect1.xlsx');
A=HNConnect (1:18,1:18);
[M,N]=size(A);
E=eye(M,N);
K=A+E;
if nargin<=3
	graph=1;
end
if graph==1
	r=.5;
	z1=(r*exp(2*pi*i*(0:6)/7)).';
	z2=(2.0*r*exp(20*pi*i*(0:8)/(9*21))).';
	z3=2.4*r*exp((pi+.1)*i/4);
	xy=zeros(M,2);
	xy([15 18 4 12 17 2 7],:)=[real(z1),imag(z1)];
	xy([9 3 1 6 16 5 14 13 10],:)=[real(z2),imag(z2)];
	xy(11,:)=[real(z3),imag(z3)];
	subplot(211);
	gplot (A,xy,'k-');
	hold on;
	plot(xy(:,1),xy(:,2),'k.','markersize',20,xy(8,1),...
	xy(8,2),'r.','markersize',20,xy(6,1),xy(6,2),...
	'b.','markersize',20);
	hold off;
	axis off;
	axis square;
	subplot(212);
	XY=randn(M,2);
	gplot (A,XY,'k-');
	hold on;
	plot(XY(:,1),XY(:,2),'k.','markersize',20,XY(8,1),...
	XY(8,2),'r.','markersize',20,XY(6,1),XY(6,2),...
	'b.','markersize',20);
	hold off;
	axis off;
	axis square;
end
x0=COVIDHist (1:18,m);
f0=find(abs(x0)<=tol);
x=COVIDHist (1:18,n);
f1=find(abs(x)<=tol);
f2=find(abs(x)>tol);
x0(f0)=0;
x(f1)=0;
T=E;
for k=f2
	T0=E;
	T0(k,:)=K(k,:).*(x(k)/x0);
	T=T0*T;
end
T0=ones(M,1);
y0=T*x0;
T0(f2)=x(f2)./y0(f2);
T=diag(T0)*T;
K=A+E;
end