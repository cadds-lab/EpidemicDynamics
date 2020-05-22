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
## function W=UACWeight(A0,A1,x0)
##
## Example:
## [x,A0]=UACPredictor(48,1e-10);
## [x,A1]=UACDescriptor(48);
## W=UACWeight(A0,A1,x(:,47:49));
## A2=[A0;A1];
## [x(:,49) A0*x(:,48) A1*x(:,48) W*A2*x(:,48)]


## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-28

function W=UACWeight(A0,A1,x0)
	[m,N]=size(x0);
	y0=x0(:,1);
	y1=x0(:,1);
	for k=1:(N-1)
		y0=[y0 A0*y0(:,k)];
		y1=[y1 A1*y1(:,k)];
	end
	for j=1:m
		c(j,:)=lsqnonneg([y0(j,:).' y1(j,:).'],x0(j,:).').';
	end
	W=[diag(c(:,1)) diag(c(:,2))];
end