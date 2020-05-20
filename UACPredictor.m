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
## function [t,x,y,A]=UACPredictor(j,tol)
##
## Example:
## [t,x,y,A]=UACPredictor(5,1e-19);

## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-28

function [t,x,y,A]=UACPredictor(j,tol,T)
	pkg load io;
	x=xlsread ('COVID19HistoryJoint.xlsx');
  	[p,m]=size(x);
	x=x(1:(p-1),:);
	[u,s,v]=svd(x,0);
  	rk=max(find(diag(s)>=tol));
	P=u(:,1:rk);	
	x0=P'*x;
	y0=x0(:,1:(m-1));
	y1=x0(:,2:m);
	for k=1:rk
		A(k,:)=lsqnonneg(y0.',y1(k,:).').';
	end
	t=1:m;
	E=@(k,n)((1:n)==k);
	K=E(j,p-1);
	y=x(:,1);
	for k=1:(T-1)
		y=[y P*A*P'*y(:,k)];
	end	
	subplot(211);
	plot(t,x(j,:),'k.-',t,y(j,1:m),'r.-');
	subplot(212);
	plot(1:T,y(j,:),'b.-');
  	A=P*A*P';
end