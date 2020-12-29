%% TPS-Warp from feature correspondences with smoothing

function [M1 M2]=TPSCoeffDiff(p,L,C,lambda,Epsilon_lambda)
if(nargin<7) %% ifndef Compute Epsilon_lambda here
    Epsilon_lambda=TPSEpsilonLambda(C,lambda);
end

[m,n]=size(p);
[np,d]=size(L);
M1=zeros(n,np);
M2=zeros(m*n,np);

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
M1(i,:)=vp';
for k=1:m
    vvec=[zeros(m+1,1)];
    vvec(k)=1;
    lpp(:,k)=[rhosd.*2.*(pvec(:,k)-C(:,k));vvec];
    M2((i-1)*m+k,:)=(Epsilon_lambda'*lpp(:,k))';
end

end

