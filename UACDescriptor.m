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
## function [Xh,T,EIHUB,EIHLB,EGHUB,EGHLB]=UACPredictor(n,r,tol)
##
## Example:
## [Xh,T,EIHUB,EIHLB,EGHUB,EGHLB]=UACPredictor(9,12,1e-12);

## Author: fredy <fredy@HPCLAB>
## Created: 2020-03-28

function [Xh,T,EIHUB,EIHLB,EGHUB,EGHLB]=UACPredictor(n,r,tol)
	pkg load io;
	Xh=xlsread ('COVID19History.xlsx');
  HNConnect=xlsread ('HNConnect0.xlsx');
	[p,m]=size(Xh);
	Xh=Xh(1:(p-1),(n+1):(r+1));
	[p,m]=size(Xh);
	Xh0=Xh(:,1:(m-1));
	[uh,sh,vh]=svd(Xh0,0);
	sh0=diag(sh);
	f=find(sh0<=tol);
	sh0(f)=tol;
	sh0=diag(sh0);
	T=Xh0\Xh(:,m);
	T=[[zeros(1,m-2);eye(m-2)] T];
	T=uh*sh0*vh'*T*(vh/sh0)*uh';
	EIHUB=Xh(:,1);
	EGHUB=Xh(:,1);
	EIHLB=EIHUB;
	EGHLB=EGHUB;	
	for k=1:(m-1)
		EIHUB = [EIHUB (Xh(:,k)>0).*ceil(T*Xh(:,k))];
		EIHLB = [EIHLB (Xh(:,k)>0).*floor(T*Xh(:,k))];
		EGHUB = [EGHUB (Xh(:,k)>0).*ceil(T*EGHUB(:,k))];
		EGHLB = [EGHLB (Xh(:,k)>0).*floor(T*EGHLB(:,k))];	
	end
end