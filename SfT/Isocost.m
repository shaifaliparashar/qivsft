function [cost,grad,Hess]=Isocost(L,p,q,C,mu,rho,Epsilon_lambda,proi,delta)
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
    [dp,~]=TPSWarpDiff(proi,delta.L,delta.C,delta.ir,delta.EpsilonLambda);
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

% Isometric constraints
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
c1=(lambda1-cdelta1).^2+(lambda2-cdelta2).^2;
c2=(lambda3-cdelta3).^2;

% Isometric constraints Term
costIso=n*rho*(sum(c1)+sum(c2));
costIsodLx=n*rho*(4*(Nx'.*(ux'*ones(1,l))')*(lambda1-cdelta1)'+4*(Ny'.*(uy'*ones(1,l))')*(lambda2-cdelta2)'+2*(Nx'.*(uy'*ones(1,l))'+Ny'.*(ux'*ones(1,l))')*(lambda3-cdelta3)');
costIsodLy=n*rho*(4*(Nx'.*(vx'*ones(1,l))')*(lambda1-cdelta1)'+4*(Ny'.*(vy'*ones(1,l))')*(lambda2-cdelta2)'+2*(Nx'.*(vy'*ones(1,l))'+Ny'.*(vx'*ones(1,l))')*(lambda3-cdelta3)');
costIsodLz=n*rho*(4*(Nx'.*(zx'*ones(1,l))')*(lambda1-cdelta1)'+4*(Ny'.*(zy'*ones(1,l))')*(lambda2-cdelta2)'+2*(Nx'.*(zy'*ones(1,l))'+Ny'.*(zx'*ones(1,l))')*(lambda3-cdelta3)');
costIsoHess=sparse(3*l,3*l);
costIsoHess(1:l,1:l)=n*rho*2*( 4*(Nx'.*(ux'*ones(1,l))')*(Nx'.*(ux'*ones(1,l))')'+ 4*(Ny'.*(uy'*ones(1,l))')*(Ny'.*(uy'*ones(1,l))')' + (Nx'.*(uy'*ones(1,l))'+Ny'.*(ux'*ones(1,l))')*(Nx'.*(uy'*ones(1,l))'+Ny'.*(ux'*ones(1,l))')');
costIsoHess(1:l,l+1:2*l)=n*rho*2*( 4*(Nx'.*(ux'*ones(1,l))')*(Nx'.*(vx'*ones(1,l))')'+ 4*(Ny'.*(uy'*ones(1,l))')*(Ny'.*(vy'*ones(1,l))')' + (Nx'.*(uy'*ones(1,l))'+Ny'.*(ux'*ones(1,l))')*(Nx'.*(vy'*ones(1,l))'+Ny'.*(vx'*ones(1,l))')');
costIsoHess(1:l,2*l+1:3*l)=n*rho*2*( 4*(Nx'.*(ux'*ones(1,l))')*(Nx'.*(zx'*ones(1,l))')'+ 4*(Ny'.*(uy'*ones(1,l))')*(Ny'.*(zy'*ones(1,l))')' + (Nx'.*(uy'*ones(1,l))'+Ny'.*(ux'*ones(1,l))')*(Nx'.*(zy'*ones(1,l))'+Ny'.*(zx'*ones(1,l))')');
costIsoHess(l+1:2*l,1:l)=costIsoHess(1:l,l+1:2*l);
costIsoHess(l+1:2*l,l+1:2*l)=n*rho*2*( 4*(Nx'.*(vx'*ones(1,l))')*(Nx'.*(vx'*ones(1,l))')'+ 4*(Ny'.*(vy'*ones(1,l))')*(Ny'.*(vy'*ones(1,l))')' + (Nx'.*(vy'*ones(1,l))'+Ny'.*(vx'*ones(1,l))')*(Nx'.*(vy'*ones(1,l))'+Ny'.*(vx'*ones(1,l))')');
costIsoHess(l+1:2*l,2*l+1:3*l)=n*rho*2*( 4*(Nx'.*(vx'*ones(1,l))')*(Nx'.*(zx'*ones(1,l))')'+ 4*(Ny'.*(vy'*ones(1,l))')*(Ny'.*(zy'*ones(1,l))')' + (Nx'.*(vy'*ones(1,l))'+Ny'.*(vx'*ones(1,l))')*(Nx'.*(zy'*ones(1,l))'+Ny'.*(zx'*ones(1,l))')');
costIsoHess(2*l+1:3*l,1:l)=costIsoHess(1:l,2*l+1:3*l);
costIsoHess(2*l+1:3*l,1+l:2*l)=costIsoHess(1+l:2*l,2*l+1:3*l);
costIsoHess(2*l+1:3*l,1+2*l:3*l)=n*rho*2*( 4*(Nx'.*(zx'*ones(1,l))')*(Nx'.*(zx'*ones(1,l))')'+ 4*(Ny'.*(zy'*ones(1,l))')*(Ny'.*(zy'*ones(1,l))')' + (Nx'.*(zy'*ones(1,l))'+Ny'.*(zx'*ones(1,l))')*(Nx'.*(zy'*ones(1,l))'+Ny'.*(zx'*ones(1,l))')');

cost=costData+costSmooth+costIso;
grad=[costDatadLx+costSmoothdLx+costIsodLx;costDatadLy+costSmoothdLy+costIsodLy;costDatadLz+costSmoothdLz+costIsodLz];
Hess=costDataHess+costSmoothHess+costIsoHess;

end
