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
## function [x,A]=UACDescriptor(ss)
##
## Example:
## A=UACDescriptor(48);

## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-28

function [x,A]=UACDescriptor(ss)
	pkg load io;
	x=xlsread ('COVID19HNHistory.xlsx');
  	K=xlsread ('HNConnect0.xlsx');
	[p,m]=size(x);
	x=x(2:p,:);
	x0=x;
	ss=min([ss m-1]);
	y0=x0(:,1:ss);
	y1=x0(:,2:(ss+1));
	E=eye(p-1);
	A=E;
	K=K+E;
	for k=1:(p-1)
		w=find(K(k,:));
		A0=E;
		A0(k,w)=lsqnonneg((y0(w,:)).',y1(k,:).').';
		A=A0*A;
	end
end