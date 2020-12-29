% TPS Warp Interpolation of points, giving derivatives
%
% SYNTAX
%   [dq,q]=TPSWarpDiff(p,L,C,lambda,EpsilonLambda)
%  
% INPUT ARGUMENTS
%   -p: a (mxn) array of mx1 vectors in the template
%   -L: TPS warp coefficients
%   -C: TPS centers in the template
%   -lambda: internal smoothing parameter. usually 1e-4
%   -EpsilonLambda: optional kernel matrix
%
% OUTPUT ARGUMENTS
%   -dq: a ((r*m)xn) vector of derivatives. The first m rows contain the
%   derivatives with respect of the first coordinate in p. Then it follows
%   the same way with the rest of coordinates
%   -q: a  rxn interpolated vector 
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

function [dq,q]=TPSWarpDiff(p,L,C,lambda,Epsilon_lambda)
if(nargin<5) %% ifndef Compute Epsilon_lambda here
    Epsilon_lambda=TPSEpsilonLambda(C,lambda);
if(nargin<4)
    lambda=1e-4
    end

end

[m,n]=size(p);
[~,d]=size(L);
q=zeros(d,n);
dq=zeros(m*d,n);
for i=1:n
% q=L*v(p);--> v(p)=epsilon_lambda*lp;
l=size(C,1);
lpp=zeros(l+m+1,m);
%Ct=[C,ones(l,1)];
% Obtaining p_vec
pvec=zeros(l,m);
for k=1:m
pvec(:,k)=p(k,i);
end
if(m>1)
distances=sum((pvec'-C').^2)';
else
    distances=((pvec'-C').^2)';
end

[rhos,rhosd]=TPSrho(distances);
lp=[rhos;p(:,i);1];
vp=Epsilon_lambda'*lp;
q(:,i)=L'*vp;

for k=1:m
    vvec=[zeros(m+1,1)];
    vvec(k)=1;
    lpp(:,k)=[rhosd.*2.*(pvec(:,k)-C(:,k));vvec];
    dq(1+d*(k-1):k*d,i)=L'*Epsilon_lambda'*lpp(:,k);
end
end





