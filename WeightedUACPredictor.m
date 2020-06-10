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
##  Example:
##  [x,A0]=UACPredictor(72,1e-10);
##  [x,A1]=UACDescriptor(72);
##  W0=UACWeight(A0,A1,x(:,70:73));
##  W1=UACWeight(A0,A1,x(:,71:73));
##  W2=UACWeight(A0,A1,x(:,72:73));
##  [t,y0]=WeightedUACPredictor(W0,A0,A1,x(:,72),16);
##  [t,y1]=WeightedUACPredictor(W1,A0,A1,x(:,72),16);
##  [t,y2]=WeightedUACPredictor(W2,A0,A1,x(:,72),16);
##  W=UACWeightedShooting(x(:,[72 73 88 89]),y0(:,[1 2 15 16]),...
##  y1(:,[1 2 15 16]),y2(:,[1 2 15 16]));
##  yt=W*[y0;y1;y2];
##  subplot(311);
##  j=6;plot([0 1 14 15],x(j,[72 73 88 89]),'c.','markersize',...
##  16,t,y0(j,:),'r.-',t,y2(j,:),'b.-',t,y1(j,:),'g.-',t,...
##  yt(j,:),'k.-');
##  axis tight
##  subplot(312);
##  j=8;plot([0 1 14 15],x(j,[72 73 88 89]),'c.','markersize',...
##  16,t,y0(j,:),'r.-',t,y2(j,:),'b.-',t,y1(j,:),'g.-',t,...
##  yt(j,:),'k.-');
##  axis tight
##  subplot(313);
##  j=19;plot([0 1 14 15],x(j,[72 73 88 89]),'c.','markersize',...
##  16,t,y0(j,:),'r.-',t,y2(j,:),'b.-',t,y1(j,:),'g.-',t,...
##  yt(j,:),'k.-');
##  axis tight

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