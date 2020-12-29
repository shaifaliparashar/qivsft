function [cost,grad,Hess]=Concost(L,p,q,C,mu,rho,Epsilon_lambda,proi,delta)
if(nargin<9)
delta=[];
end

l=size(C,1);
n=size(p,2);
d=size(q,1);
% for each p=(x y) obtain the v(p)
PSI=q';
N=zeros(n,l);
Nx=zeros(n,l);
Ny=zeros(n,l);
for i=1:n
    pp=p(:,i);
    % Obtaining p_vec
    pvec=zeros(l,2);
    pvec(:,1)=pp(1);
    pvec(:,2)=pp(2);
    distances=sum((pvec'-C').^2)';
    [rhos,rhosd]=TPSrho(distances);
    lp=[rhos;pp;1];
    N(i,:)=Epsilon_lambda'*lp;
    lpx=[rhosd.*2.*(pvec(:,1)-C(:,1));1;0;0];
    lpy=[rhosd.*2.*(pvec(:,2)-C(:,2));0;1;0];
    Nx(i,:)=Epsilon_lambda'*lpx;
    Ny(i,:)=Epsilon_lambda'*lpy;
end
Epsilon_bar=Epsilon_lambda(1:size(Epsilon_lambda,1)-3,:);%inv(Klambda)*(eye(l)-Ct*A);
ZTZ=8.*pi.*Epsilon_bar;
% Data Term
spdq1=sparse(diag(q(1,:)));
spdq2=sparse(diag(q(2,:)));
spdq2s=sparse(diag(q(2,:).^2));
spdq1s=sparse(diag(q(1,:).^2));
costData=(N*L(:,1)-spdq1*N*L(:,3))'*(N*L(:,1)-spdq1*N*L(:,3))+(N*L(:,2)-spdq2*N*L(:,3))'*(N*L(:,2)-spdq2*N*L(:,3));
costDatadLx=2*N'*N*L(:,1)-2*N'*spdq1*N*L(:,3);
costDatadLy=2*N'*N*L(:,2)-2*N'*spdq2*N*L(:,3);
costDatadLz=2*N'*spdq1s*N*L(:,3)-2*N'*spdq1*N*L(:,1)+2*N'*spdq2s*N*L(:,3)-2*N'*spdq2*N*L(:,2);
costDataHess=sparse(3*l,3*l);
costDataHess(1:l,1:l)=2*N'*N;
costDataHess(1:l,2*l+1:3*l)=-2*N'*spdq1*N;
costDataHess(l+1:2*l,l+1:2*l)=2*N'*N;
costDataHess(l+1:2*l,2*l+1:3*l)=-2*N'*spdq2*N;
costDataHess(2*l+1:3*l,1:l)=-2*N'*spdq1*N;
costDataHess(2*l+1:3*l,l+1:2*l)=-2*N'*spdq2*N;
costDataHess(2*l+1:3*l,2*l+1:3*l)=2*N'*spdq1s*N+2*N'*spdq2s*N;

% Smooth Term
costSmooth=trace(n*mu*L'*ZTZ*L);
costSmoothdLx=2*n*mu*ZTZ*L(:,1);
costSmoothdLy=2*n*mu*ZTZ*L(:,2);
costSmoothdLz=2*n*mu*ZTZ*L(:,3);
costSmoothHess=2*n*mu*blkdiag(ZTZ,ZTZ,ZTZ);

nroi=size(proi,2);
 if(~isempty(delta))
    [dp,~]=TPSWarpDiff(proi,options.delta.L,options.delta.C,options.delta.ir,options.delta.EpsilonLambda);
    cdelta1=dp(1,:).^2+dp(2,:).^2+dp(3,:).^2;
    cdelta2=dp(4,:).^2+dp(5,:).^2+dp(6,:).^2;
    cdelta3=dp(1,:).*dp(4,:)+dp(2,:).*dp(5,:)+dp(3,:).*dp(6,:);
 else
     cdelta1=1;cdelta2=1;cdelta3=0;
 end            


for i=1:nroi
    pp=proi(:,i);
    % Obtaining p_vec
    pvec=zeros(l,2);
    pvec(:,1)=pp(1);
    pvec(:,2)=pp(2);
    distances=sum((pvec'-C').^2)';
    [rhos,rhosd]=TPSrho(distances);
    lp=[rhos;pp;1];
    N(i,:)=Epsilon_lambda'*lp;
    lpx=[rhosd.*2.*(pvec(:,1)-C(:,1));1;0;0];
    lpy=[rhosd.*2.*(pvec(:,2)-C(:,2));0;1;0];
    Nx(i,:)=Epsilon_lambda'*lpx;
    Ny(i,:)=Epsilon_lambda'*lpy;
end

% Conformal constraints
c1=0;c2=0;
ux=L(:,1)'*Nx';
uy=L(:,1)'*Ny';
vy=L(:,2)'*Ny';
vx=L(:,2)'*Nx';
zy=L(:,3)'*Ny';
zx=L(:,3)'*Nx'; 
lambda1=ux.*ux+vx.*vx+zx.*zx;
lambda2=uy.*uy+vy.*vy+zy.*zy;
lambda3=ux.*uy+vx.*vy+zx.*zy;
c1=(cdelta2.*lambda1-cdelta1.*lambda2).^2;
c2=(cdelta3.*lambda1-cdelta1.*lambda3).^2;

% Isometric constraints Term
costIso=n*rho*(sum(c1)+sum(c2));
vecx1=2.*(Nx'.*((cdelta2.*ux)'*ones(1,l))'-Ny'.*((cdelta1.*uy)'*ones(1,l))');
vecx2=(2.*Nx'.*((cdelta3.*ux)'*ones(1,l))'-Nx'.*((cdelta1.*uy)'*ones(1,l))'-Ny'.*((cdelta1.*ux)'*ones(1,l))');
costIsodLx=n*rho*(2*vecx1*(cdelta2.*lambda1-cdelta1.*lambda2)'+2*vecx2*(cdelta3.*lambda1-cdelta1.*lambda3)');
vecy1=2.*(Nx'.*((cdelta2.*vx)'*ones(1,l))'-Ny'.*((cdelta1.*vy)'*ones(1,l))');
vecy2=(2.*Nx'.*((cdelta3.*vx)'*ones(1,l))'-Nx'.*((cdelta1.*vy)'*ones(1,l))'-Ny'.*((cdelta1.*vx)'*ones(1,l))');
costIsodLy=n*rho*(2*vecy1*(cdelta2.*lambda1-cdelta1.*lambda2)'+2*vecy2*(cdelta3.*lambda1-cdelta1.*lambda3)');
vecz1=2.*(Nx'.*((cdelta2.*zx)'*ones(1,l))'-Ny'.*((cdelta1.*zy)'*ones(1,l))');
vecz2=(2.*Nx'.*((cdelta3.*zx)'*ones(1,l))'-Nx'.*((cdelta1.*zy)'*ones(1,l))'-Ny'.*((cdelta1.*zx)'*ones(1,l))');
costIsodLz=n*rho*(2*vecz1*(cdelta2.*lambda1-cdelta1.*lambda2)'+2*vecz2*(cdelta3.*lambda1-cdelta1.*lambda3)');

costIsoHess=sparse(3*l,3*l);
costIsoHess(1:l,1:l)=n*rho*2*( vecx1*vecx1'+ vecx2*vecx2');
costIsoHess(1:l,l+1:2*l)=n*rho*2*(vecx1*vecy1'+vecx2*vecy2');
costIsoHess(1:l,2*l+1:3*l)=n*rho*2*(vecx1*vecz1'+vecx2*vecz2');
costIsoHess(l+1:2*l,1:l)=costIsoHess(1:l,l+1:2*l);
costIsoHess(l+1:2*l,l+1:2*l)=n*rho*2*(vecy1*vecy1'+vecy2*vecy2');
costIsoHess(l+1:2*l,2*l+1:3*l)=n*rho*2*(vecy1*vecz1'+vecy2*vecz2');
costIsoHess(2*l+1:3*l,1:l)=costIsoHess(1:l,2*l+1:3*l);
costIsoHess(2*l+1:3*l,1+l:2*l)=costIsoHess(1+l:2*l,2*l+1:3*l);
costIsoHess(2*l+1:3*l,1+2*l:3*l)=n*rho*2*(vecz1*vecz1'+vecz2*vecz2');

cost=costData+costSmooth+costIso;
grad=[costDatadLx+costSmoothdLx+costIsodLx;costDatadLy+costSmoothdLy+costIsodLy;costDatadLz+costSmoothdLz+costIsodLz];
Hess=costDataHess+costSmoothHess+costIsoHess;

end
