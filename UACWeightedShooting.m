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
## function W=UACWeightedShooting(x,yl,ym,yh)
##

## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-28

function W=UACWeightedShooting(x,yl,ym,yh)
	[m,N]=size(x);
	for j=1:m
		c(j,:)=lsqnonneg([yl(j,:).' ym(j,:).' yh(j,:).'],x(j,:).').';
	end
	W=[diag(c(:,1)) diag(c(:,2)) diag(c(:,3))];
end