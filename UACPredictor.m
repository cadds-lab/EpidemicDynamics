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
## function [x,A]=UACPredictor(ss,tol)
##
## Example:
## A=UACPredictor(48,1e-19);

## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-28

function [x,A]=UACPredictor(ss,tol)
	pkg load io;
	%x=xlsread ('COVID19HNHistory.xlsx');
	x=xlsread ('COVID19HNHistoryFull.xlsx');
  	[p,m]=size(x);	
	ss=min([ss m-1]);
	x=x(2:p,1:(ss+1));
	[u,s,v]=svd(x,0);
  	rk=max(find(diag(s)>=tol));
	P=u(:,1:rk);	
	x0=P'*x;
	y0=x0(:,1:ss);
	y1=x0(:,2:(ss+1));
	for k=1:rk
		A(k,:)=lsqnonneg(y0.',y1(k,:).').';
	end
  	A=P*A*P';
end