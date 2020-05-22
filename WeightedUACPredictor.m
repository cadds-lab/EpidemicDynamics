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
## function [t,x]=WeightedUACPredictor(W,A0,A1,x0,T)
##
## Example:
## [x,A0]=UACPredictor(48,1e-10);
## [x,A1]=UACDescriptor(48);
## W0=UACWeight(A0,A1,x(:,46:49));
## W1=UACWeight(A0,A1,x(:,47:49));
## W2=UACWeight(A0,A1,x(:,48:49));
## [t,y0]=WeightedUACPredictor(W0,A0,A1,x(:,48),25);
## [t,y1]=WeightedUACPredictor(W1,A0,A1,x(:,48),25);
## [t,y2]=WeightedUACPredictor(W2,A0,A1,x(:,48),25);
## W=UACWeightedShooting(x(:,48:51),y0(:,[1 2 24 25]),...
## y1(:,[1 2 24 25]),y2(:,[1 2 24 25]));
## yt=W*[y0;y1;y2];
## subplot(311);
## j=6;plot([0 1 23 24],x(j,[48 49 50 51]),'c.','markersize',...
## 16,t,y0(j,:),'r.-',t,y2(j,:),'b.-',t,y1(j,:),'g.-',t,...
## yt(j,:),'k.-');
## subplot(312);
## j=8;plot([0 1 23 24],x(j,[48 49 50 51]),'c.','markersize',...
## 16,t,y0(j,:),'r.-',t,y2(j,:),'b.-',t,y1(j,:),'g.-',t,...
## yt(j,:),'k.-');
## subplot(313);
## j=19;plot([0 1 23 24],x(j,[48 49 50 51]),'c.','markersize',...
## 16,t,y0(j,:),'r.-',t,y2(j,:),'b.-',t,y1(j,:),'g.-',t,...
## yt(j,:),'k.-');



## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-28

function [t,x]=WeightedUACPredictor(W,A0,A1,x0,T)
	x=x0(:,1);
	A=[A0;A1];
	for k=1:(T-1)
		x=[x W*A*x(:,k)];
	end
	t=0:(T-1);
end