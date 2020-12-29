% TPS Coefficients from feature correspondences. It uses Linear Least Squares.
%
% SYNTAX
%   L=TPSWfromfeatures(x1,x2,C,mu,lambda,EpsilonLambda)
%  
% INPUT ARGUMENTS
%   -x1: a (mxn) array of mx1 vectors in the template
%   -x2: a (rxn) array of rx1 vectors as output of the warp
%   -C: TPS centers in the template
%   -mu: smoothing parameter
%    -lambda: internal smoothing parameter. usually 1e-4
%   -EpsilonLambda: optional kernel matrix
%
% OUTPUT ARGUMENTS
%   -L: warp coefficients
%
% (c) 2013, Daniel Pizarro. dani.pizarro@gmail.com 

% TPS is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% TPS is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


function [L,ZTZ]=TPSWfromfeatures(x1,x2,C,mu,lambda,Epsilon_lambda)
if(nargin<6) %% ifndef Compute Epsilon_lambda here
    Epsilon_lambda=TPSEpsilonLambda(C,lambda);
end
[l,m]=size(C);
% for each p=(x y) obtain the v(p)
PSI=[];
N=[];

for i=[1:size(x1,2)]
    p=x1(:,i);
    q=x2(:,i);
% q=L*v(p);--> v(p)=epsilon_lambda*lp;

% Obtaining p_vec
pvec=zeros(l,m);
for k=1:m
pvec(:,k)=p(k);
end
if(m>1)
distances=sum((pvec'-C').^2)';
else
    distances=((pvec'-C').^2)';
end
rhos=TPSrho(distances);
lp=[rhos;p;1];
N=[N,Epsilon_lambda'*lp];
PSI=[PSI,q];
end
N=N';
PSI=PSI';   
mf=length(x1);
Epsilon_bar=Epsilon_lambda(1:size(Epsilon_lambda,1)-(m+1),:);%inv(Klambda)*(eye(l)-Ct*A);
ZTZ=8.*pi.*Epsilon_bar;
L=(((N'*N)+mf*mu*ZTZ))\(N'*PSI);
