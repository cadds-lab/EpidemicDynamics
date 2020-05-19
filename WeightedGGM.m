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
## function [Y,t]=WeightedGGM(department,graph)
##
## Example: [Y,t]=WeightedGGM(19);
## 

## Author: fredy <fredy@HPCLAB>
## Created: 2020-04-10

function [Y,t]=WeightedGGM(department,graph)

if nargin==0, department=19;end
if nargin<=1, graph=1;end

pkg load io;
Xh = xlsread ('COVID19HistoryJoint.xlsx');
X=Xh(department,:);

if sum(X)==0 
	Y1=Y2=Y3=Y=X;
	t=1:length(X);
else
	[Y1,t,D,ts]=SDiffLB(X,1);
	[Y2,t]=EMod(X);
	[Y3,t,D,ts]=SDiffUB(X,1);
	Y=[Y1;Y2;Y3];
	w=Y'\X';
	Y=w'*Y;
end
if graph==1
	plot(t,X,'k.-','markersize',...
	12,t,Y1,...
	'b.-','markersize',12,t,Y2,...
	'c.-','markersize',12,t,Y3,...
	'r.-','markersize',12,t,Y,...
	'g.-','markersize',12);
	legend('C(t)','C_0(t)',...
	'C_1(t)','C_2(t)','C_{w}(t)');
	grid on;
end
end

function [Yt,t,D,ts]=SDiffLB(S,s)
Ls=length(S);
t=0:(Ls-1);
ds=diff(S);
ff=find(ds>0);
D=ds(ff(1));
ts=t(ff(1));
Lff=length(ff);
for k=2:Lff-1
	if s==0,if sum(ds(ff(k))<=...
	ds(ff((k+1):Lff)))>=Lff-k-1 ... 
	&& ds(ff(k))<=ds(ff(k+1)), ...
	D=[D ds(ff(k))];ts=[ts ...
	t(ff(k))];end;end
	if s==1,if sum(ds(ff(k))<...
	ds(ff((k+1):Lff)))>=Lff-k-1 ...
	&& ds(ff(k))<ds(ff(k+1)), D=[D ...
	ds(ff(k))];ts=[ts t(ff(k))];end;
	end
end
D=[D ds(ff(Lff))];
ts=[ts t(ff(Lff))];
LD=log(D);
p=polyfit(ts,LD,1);
Yt=S(1);
for k=1:(Ls-1) 
	Yt=[Yt Yt(k)+...
	exp(polyval(p,k-1))];
end
end

function [Yt,t,D,ts]=SDiffUB(S,s)
Ls=length(S);
t=0:(Ls-1);
ds=diff(S);
ff=find(ds>0);
D=ds(ff(1));
ts=t(ff(1));
for k=2:length(ff)
	if s==0,if sum(ds(ff(k))>=...
	ds(ff(1:(k-1))))>=k-1, D=...
	[D ds(ff(k))];ts=[ts ...
	t(ff(k))];end;end
	if s==1,if sum(ds(ff(k))>...
	ds(ff(1:(k-1))))>=k-1, D=...
	[D ds(ff(k))];ts=[ts ...
	t(ff(k))];end;end
end
LD=log(D);
p=polyfit(ts,LD,1);
Yt=S(1);
for k=1:(Ls-1) 
	Yt=[Yt Yt(k)+...
	exp(polyval(p,k-1))];
end
end

function [Yt,t]=EMod(S)
Ls=length(S);
t=0:(Ls-1);
D=diff(S);
ff=find(D>0);
ts=t(ff);
LD=log(D(ff));
p=polyfit(ts,LD,1);
Yt=S(1);
for k=1:(Ls-1) 
	Yt=[Yt Yt(k)+...
	exp(polyval(p,k-1))];
end
end