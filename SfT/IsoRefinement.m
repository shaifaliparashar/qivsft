function phi=IsoRefinement(p,q,phi,options)
if(nargin<3)
    disp('Error: 2xn arrays p and q and initialization phi are needed...');
    coeff=[];
    out=[];
    return
elseif(nargin==3)
    [options,error]=ProcessArgs(p,q,phi);
else
    [options,error]=ProcessArgs(p,q,phi,options);
end


%% L-M refinement

iter=0; nu=0.01; found=0;delta=options.tol;
X=phi.L;
[cost,costgrad,Hess]=Isocost(X,p,q,phi.C,phi.er,options.isoer,phi.EpsilonLambda,phi.p,options.delta);

while(found==0 & iter<options.maxiter)

    costant=cost;
    costgradant=costgrad;
    Hessant=Hess;
    Xant=X;
    while(rcond(Hess+nu.*eye(size(Hess)))<1e-12)
    nu=nu*2;
    end
    XDelta=reshape(-inv(Hess+nu.*eye(size(Hess)))*costgrad,size(X));
[cost,costgrad,Hess]=Isocost(X+XDelta,p,q,phi.C,phi.er,options.isoer,phi.EpsilonLambda,phi.p,options.delta);
    if(cost<costant)
   %    disp(sprintf('cost_ant=%f',costant))
        X=X+XDelta;
        if(nu>0.0000001)
        nu=nu/2;     
        end
        found=(norm(X(:)-Xant(:))<delta);
        %found=norm(costgrad)<delta;
    else
        nu=nu*10;
        cost=costant;
        costgrad=costgradant;
        Hess=Hessant;
    end
found=(found|(nu>1e10));
if(options.verbose)
 disp(sprintf('[REISO-PHI]iter=%d;cost=%f;nu=%f;diff=%f',iter,cost,nu,norm(costgrad)))
end
    iter=iter+1;
end
phi.L=X;
end



function [options,error]=ProcessArgs(p,q,phi,options)
% hard coded defaults
nC=10;ir=1e-4;er=0.55;
%
error=[];
[d,n]=size(p);
[d2,n2]=size(q);
if(d<2 || d2<2 || n2~=n || n<3)
    error='Point arrays with mismatched dimmensions or few points given...';
    options=[];
    return
end
if(nargin<3)
options=[];
end
if(~isfield(options,'maxiter'))
        options.maxiter=40;
end
if(~isfield(options,'isoer'))
        options.isoer=10;
end
if(~isfield(options,'tol'))
        options.tol=1e-4;
end
if(~isfield(options,'verbose'))
        options.verbose=1;
end
if(~isfield(options,'delta'))
        options.delta=[];
end
end


