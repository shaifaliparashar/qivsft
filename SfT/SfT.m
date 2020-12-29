
% Shape from Template
%
% SYNTAX
%   [out]=Sft(p,q,options)
%
% INPUT ARGUMENTS
%   -p: (2xn) array of n 2D points in the flat template
%   -q: (2xn) array of b 2D points in the image. q must be normalised with
%   the intrinsic parameter matrix of the camera.
%   -options: structure with the following fields
%           'eta': struct with eta (template to image warp) parameters
%                  eta.ir: internal smoothing. default=1e-4;
%                  eta.er: external smoothing. default=0.55;
%                  eta.nC: nC^2 control centers
%           'phi': struct with phi (template to shape warp) parameters
%                   phi.ir: internal smoothing, default=1e-4;
%                   phi.er: external smoothing, default=0.55;
%                   phi.nC: nC^2 control centers
%           'verbose': 1 --> gives debug information
%           'KLims': Rectangle bounds of the template.
%                     KLims=[umin,umax,vmin,vmax]
%           'method': method used for shape estimation
%                    'AnIso' --> Analytical solution for Isometric
%                    Deformations.
%                    'ReIso' --> Analytical solution + refinement for
%                    Isometric Deformations
%                    'AnCon' --> Analytical solution for Conformal
%                    Deformations
%                    'ReCon' --> Analytical solution + refinements for
%                    Conformal Deformations
%           'NGridx': number of grid points in x used to sample the template (default=50)
%           'NGridy': number of grid points in y used to sample the template (default=50).
%           'maxiter': number of max iterations in case of refinement
%           (default=40)
%           'delta': (Warning: delta must be provided only if the template is 3D !)
%            warp from flat template to 3D template. Structure with TPS parameters
%                  delta.ir: internal smoothing
%                  delta.er: external smoothing
%                  delta.nC: number of control centers
%                  delta.C: control centers
%                  delta.L: warp coefficients
%                  delta.EpsilonLambda: TPS kernel matrix
%           'phigth': (Warning: Only for choosing closest solution after AnCon or ReCon !) ground truth
%            phi warp from flat template to 3D shape. Structure with TPS parameters
%                  phigth.ir: internal smoothing
%                  phigth.er: external smoothing
%                  phigth.nC: number of control centers
%                  phigth.C: control centers
%                  phigth.L: warp coefficients
%                  phigth.EpsilonLambda: TPS kernel matrix
%
% OUTPUT ARGUMENTS
%       out: output structure with the following fields:
%           'phi': solution to shape. structure with TPS parameters
%                  phi.ir: internal smoothing
%                  phi.er: external smoothing
%                  phi.nC: number of control centers
%                  phi.C: control centers
%                  phi.L: warp coefficients
%                  phi.EpsilonLambda: TPS kernel matrix
%                  phi.p: grid of points in the template
%                  phi.Q: grid of points in shape space.
%       NOTE: if method='ReCon' or 'AnCon' and we don't provide groundtruth
%       warp phi is a cell array with all solutions found.
%           'eta': registration warp. structure with TPS parameters
%                  eta.ir: internal smoothing
%                  eta.er: external smoothing
%                  eta.nC: number of control centers
%                  eta.C: control centers
%                  eta.L: warp coefficients
%                  eta.EpsilonLambda: TPS kernel matrix
%                  eta.p: grid of points in the template
%                  eta.q: grid of points in the image.
%
%   This code is partially based on the work of
%   [Bartoli et. al 2012]On Template-Based Reconstruction from a Single View:
%   Analytical Solutions and Proofs of Well-Posedness for
%   Developable, Isometric and Conformal Surfaces
%
%   (c) 2013, Adrien Bartoli and Daniel Pizarro. dani.pizarro@gmail.com
%
% Sft is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Sft is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

% SfT Shape from Template source code
function [out]=SfT(p,q,options)
if(nargin<2)
    disp('Error: 2xn arrays p1 and p2 are needed...');
    coeff=[];
    out=[];
    return
elseif(nargin==2)
    [options,error]=ProcessArgs(p,q);
else
    [options,error]=ProcessArgs(p,q,options);
end


% Get Warp Centers
C=TPSGenerateCenters(options.eta.nC,options.KLims+1e-3.*[-1,1,-1,1]);
% Precompute EpsilonLambda matrix
EpsilonLambda=TPSEpsilonLambda(C,options.eta.ir);
% Get warp parameters from features
L=TPSWfromfeatures(p,q,C,options.eta.er,options.eta.ir,EpsilonLambda);
% Warp points p and get warp reprojection error
[~,qw]=TPSWarpDiff(p,L,C,options.eta.ir,EpsilonLambda);
error=sqrt(mean((qw(1,:)-q(1,:)).^2+(qw(2,:)-q(2,:)).^2));
if(options.verbose)
    %Visualize Point Registration Error
    [xv,yv]=meshgrid(linspace(options.KLims(1),options.KLims(2),20),linspace(options.KLims(3),options.KLims(4),20));
    [~,qv]=TPSWarpDiff([xv(:)';yv(:)'],L,C,options.eta.ir,EpsilonLambda);
    disp([sprintf('[ETA] Internal Rep error = %f',error)]);
    figure;
    plot(q(1,:),q(2,:),'ro');
    hold on;
    plot(qw(1,:),qw(2,:),'b*');
    mesh(reshape(qv(1,:),size(xv)),reshape(qv(2,:),size(xv)),zeros(size(xv)));
    hold off;
    figure;
    plot(p(1,:),p(2,:),'r*');
    hold on;
    plot(C(:,1),C(:,2),'bo');
end
out.eta.p=p;
out.eta.q=q;
out.eta.L=L;
out.eta.C=C;
out.eta.EpsilonLambda=EpsilonLambda;
out.eta.ir=options.eta.ir;
out.eta.er=options.eta.er;

switch(options.method)
    case {'AnIso','ReIso'}% Isometric Solution
        if(~(isfield(options.phi,'L') & strcmp(lower(options.method),'reiso'))) % Is user giving an initialization of phi ?
            % Create a grid of points to go to 3D
            NRoi=options.NGridy*options.NGridx;
            [xroi,yroi]=meshgrid(linspace(options.KLims(1),options.KLims(2),options.NGridx),linspace(options.KLims(3),options.KLims(4),options.NGridy));
            proi=[xroi(:)';yroi(:)'];
            % Get Derivatives
            [dq,qw]=TPSWarpDiff(proi,L,C,options.eta.ir,EpsilonLambda);
            % Get phi points
            gamma=zeros(1,NRoi);
            Q=zeros(3,NRoi);
            if(isfield(options,'delta'))
                [dp,~]=TPSWarpDiff(proi,options.delta.L,options.delta.C,options.delta.ir,options.delta.EpsilonLambda);
            end
            for i=1:NRoi
                Jdelta=eye(2);
                if(isfield(options,'delta'))
                    Jdelta=[dp(1:3,i)';dp(4:6,i)']';
                end
                eta=[qw(1,i);qw(2,i)];
                Jeta=[dq(1,i),dq(3,i);dq(2,i),dq(4,i)];
                M=(Jdelta'*Jdelta)/(Jeta'*(eye(2)-(eta*eta')./(eta'*eta+1))*Jeta);
                eigM=eig(M);
                gamma(i)=(sqrt(min(eigM)));
                Q(1,i)=gamma(i)*eta(1);
                Q(2,i)=gamma(i)*eta(2);
                Q(3,i)=gamma(i);
            end
            
            % Get Warp Centers
            C=TPSGenerateCenters(options.phi.nC,options.KLims+1e-3*[-1,1,-1,1]);
            % Precompute EpsilonLambda matrix
            EpsilonLambda=TPSEpsilonLambda(C,options.phi.nC);
            L2=TPSWfromfeatures(proi,Q,C,options.phi.er,options.phi.ir,EpsilonLambda);
            [~,Qw]=TPSWarpDiff(proi,L2,C,options.phi.ir,EpsilonLambda);
            error=sqrt(mean((Qw(1,:)-Q(1,:)).^2+(Qw(2,:)-Q(2,:)).^2+(Qw(3,:)-Q(3,:)).^2));
            if(options.verbose)
                %Visualize Point Registration Error
                disp([sprintf('[PHI] Internal Rep error = %f',error)]);
            end
            out.phi.Q=Q;
            out.phi.p=proi;
            out.phi.L=L2;
            out.phi.C=C;
            out.phi.EpsilonLambda=EpsilonLambda;
            out.phi.ir=options.phi.ir;
            out.phi.er=options.phi.er;
        else
            out.phi=options.phi;
        end
        
        switch(options.method)
            case 'ReIso'
                out.init.phi=options.phi;
                phi=IsoRefinement(p,q,out.phi,options);
                out.phi=phi;
        end
        
    case {'AnCon','ReCon'}
        if(~(isfield(options.phi,'L') & strcmp(lower(options.method),'recon'))) % Is user providing with an initialization of phi ?
            % Create a grid of points to go to 3D
            NRoi=options.NGridy*options.NGridx;
            [xroi,yroi]=meshgrid(linspace(options.KLims(1),options.KLims(2),options.NGridx),linspace(options.KLims(3),options.KLims(4),options.NGridy));
            proi=[xroi(:)';yroi(:)'];
            % Get Derivatives
            [dq,qw]=TPSWarpDiff(proi,L,C,options.eta.ir,EpsilonLambda);
            if(isfield(options,'phigth'))
                [dgth,~]=TPSWarpDiff(proi,options.phigth.L,options.phigth.C,options.phigth.ir,options.phigth.EpsilonLambda);
            else
                dgth=[];
            end
            
            
            % Get phi points
            gamma=zeros(1,NRoi);
            Jmu=zeros(2,NRoi);
            Q=zeros(3,NRoi);
            I=zeros(options.NGridy,options.NGridx);
            IU=I;IV=I;Ig=I;
            if(isfield(options,'delta'))
                [dp,~]=TPSWarpDiff(proi,options.delta.L,options.delta.C,options.delta.ir,options.delta.EpsilonLambda);
            end
            for i=1:NRoi
                Jdelta=eye(2);
                if(isfield(options,'delta'))
                    Jdelta=[dp(1:3,i)';dp(4:6,i)']';
                end
                eta=[qw(1,i);qw(2,i)];
                Jeta=[dq(1,i),dq(3,i);dq(2,i),dq(4,i)]  ;
                Vt=chol((Jdelta'*Jdelta));V=Vt';
                A=(inv(V)*((Jeta'*(eye(2)-(eta*eta')./(eta'*eta+1))*Jeta))*inv(V'))./(eta'*eta+1);
                [eigvA,eigA,~]=svd(A);
                eiglist=diag(eigA);
                eigvA=eigvA*sign(eigvA(1,2));
                Jmu(:,i)=(sqrt(abs(eigA(1,1)-eigA(2,2))))*V*eigvA(:,2);
                if(length(dgth)>0)
                    if(sign(dgth(6,i))~=sign(Jmu(2,i)))
                        Ig(i)=1;
                    end
                end
                I(i)=abs(Jmu(1,i))+abs(Jmu(2,i));
            end
            % Get analytical derivatives of |Jmu| using a TPS
            C=TPSGenerateCenters(options.phi.nC,options.KLims+1e-3*[-1,1,-1,1]);
            EpsilonLambda=TPSEpsilonLambda(C,options.phi.nC);
            Ld=TPSWfromfeatures(proi,I(:)',C,options.phi.er*10,options.phi.ir,EpsilonLambda);
            [LdD,~]=TPSWarpDiff(proi,Ld,C,options.phi.ir,EpsilonLambda);
            IU(:)=LdD(1,:);
            IV(:)=LdD(2,:);
            [Bwr,nreig,ncombinations]=getcombinations(IU,0);
            [Bwr2,nreig2,ncombinations2]=getcombinations(IV,0);
            kc3=1;
            scoremax=0;
            for kc2=1:size(ncombinations2,1)
                mask2=ones(size(Bwr2));
                for kr=1:nreig2
                    if(ncombinations2(kc2,kr)=='1')
                        mask2(Bwr2==(kr))=-1;
                    end
                end
                
                for kc=1:size(ncombinations,1)
                    mask=ones(size(Bwr));
                    for kr=1:nreig
                        if(ncombinations(kc,kr)=='1')
                            mask(Bwr==(kr))=-1;
                        end
                    end
                    if(length(dgth)==0)
                        Jmu2=Jmu;
                        Jmu2(1,:)=(mask(:)').*(mask2(:)').*Jmu(1,:);
                        Jmu2(2,:)=(mask(:)').*(mask2(:)').*Jmu(2,:);
                        Lmu=TPSIntegration(proi,Jmu2,C,options.phi.er/100,options.phi.ir,EpsilonLambda);
                        [~,muw]=TPSWarpDiff(proi,Lmu,C,options.phi.ir,EpsilonLambda);
                        gamma=exp(muw)./sqrt(1+qw(1,:).^2+qw(2,:).^2);
                        Q(1,:)=gamma.*qw(1,:);
                        Q(2,:)=gamma.*qw(2,:);
                        Q(3,:)=gamma;
                        % Get Warp Centers
                        L2=TPSWfromfeatures(proi,Q,C,options.phi.er,options.phi.ir,EpsilonLambda);
                        %[~,Qw]=TPSWarpDiff(proi,L2,C,options.phi.ir,EpsilonLambda);
                        out.phi{kc3}.Q=Q;
                        out.phi{kc3}.p=proi;
                        out.phi{kc3}.L=L2;
                        out.phi{kc3}.C=C;
                        out.phi{kc3}.EpsilonLambda=EpsilonLambda;
                        out.phi{kc3}.ir=options.phi.ir;
                        out.phi{kc3}.er=options.phi.er;
                        out.phi{kc3}.I=I;
                        switch(options.method)
                            case 'ReCon'
                                out.init.phi{kc3}=out.phi{kc3};
                                phi=ConRefinement(p,q,out.phi{kc3},options);
                                out.phi{kc3}=phi;
                        end
                        kc3=kc3+1;
                    else
                        Ir=I;
                        Ir(:)=((mask(:)).*(mask2(:)));
                        
                        scorei=sum(Ir(:).*(-Ig(:).*2+1));
                        if(scorei>scoremax)
                            scoremax=scorei;
                            maskt=(mask(:)').*(mask2(:)');
                        end
                        
                    end
                    
                end
            end
            if(length(dgth)>0)
                kc3=1;
                Jmu2=Jmu;
                Jmu2(1,:)=maskt.*Jmu(1,:);
                Jmu2(2,:)=maskt.*Jmu(2,:);
                Lmu=TPSIntegration(proi,Jmu2,C,options.phi.er/100,options.phi.ir,EpsilonLambda);% /100
                [~,muw]=TPSWarpDiff(proi,Lmu,C,options.phi.ir,EpsilonLambda);
                gamma=exp(muw)./sqrt(1+qw(1,:).^2+qw(2,:).^2);
                Q(1,:)=gamma.*qw(1,:);
                Q(2,:)=gamma.*qw(2,:);
                Q(3,:)=gamma;
                % Get Warp Centers
                L2=TPSWfromfeatures(proi,Q,C,options.phi.er,options.phi.ir,EpsilonLambda);
                %[~,Qw]=TPSWarpDiff(proi,L2,C,options.phi.ir,EpsilonLambda);
                out.phi.Q=Q;
                out.phi.p=proi;
                out.phi.L=L2;
                out.phi.C=C;
                out.phi.EpsilonLambda=EpsilonLambda;
                out.phi.ir=options.phi.ir;
                out.phi.er=options.phi.er;
                out.phi.I=I;
                switch(options.method)
                    case 'ReCon'
                        out.init.phi{kc3}=out.phi{kc3}
                        phi=ConRefinement(p,q,out.phi{kc3},options);
                        out.phi{kc3}=phi;
                end
            end
        else
            kc3=1;
            phi=options.phi;
            switch(options.method)
                case 'ReCon'
                    out.init.phi=options.phi;
                    phi=ConRefinement(p,q,options.phi,options);
                    out.phi=phi;
            end
        end
        
        
    otherwise
        
        disp('Error: method not recognized');
        
end
%




end



function [options,error]=ProcessArgs(p1,p2,options)
% hard coded defaults
nC=10;ir=1e-3;er=0.55;
%
error=[];
[d,n]=size(p1);
[d2,n2]=size(p2);
if(d<2 || d2<2 || n2~=n || n<3)
    error='Point arrays with mismatched dimmensions or few points given...';
    options=[];
    return
end
if(nargin<3)
    options=[];
end
if(~isfield(options,'eta'))
    options.eta.er=er;
    options.eta.ir=ir;
    options.eta.nC=nC;
else
    if(~isfield(options.eta,'er'))
        options.eta.er=er;
    end
    if(~isfield(options.eta,'ir'))
        options.eta.ir=ir;
    end
    if(~isfield(options.eta,'nC'))
        options.eta.nC=nC;
    end
end

if(~isfield(options,'phi'))
    options.phi.er=er;
    options.phi.ir=ir;
    options.phi.nC=nC;
else
    if(~isfield(options.phi,'er'))
        options.phi.er=er;
    end
    if(~isfield(options.phi,'ir'))
        options.phi.ir=ir;
    end
    if(~isfield(options.phi,'nC'))
        options.phi.nC=nC;
    end
end



if(~isfield(options,'verbose'))
    options.verbose=0;
end
if(~isfield(options,'KLims'))
    umin=min(p1(1,:));
    umax=max(p1(1,:));
    vmin=min(p1(2,:));
    vmax=max(p1(2,:));
    options.KLims=[umin,umax,vmin,vmax];
end
if(~isfield(options,'method'))
    options.method='AnIso';
end
if(~isfield(options,'NGridx'))
    options.NGridx=50;
end
if(~isfield(options,'NGridy'))
    options.NGridy=50;
end

end
