function [Q,J_delta,J_phi]=create_template(p_roi,P_roi,q,Q_roi, e,p,scale,KK)

delta.nC=9; % control centers of warp delta
umin=min(p_roi(1,:));umax=max(p_roi(1,:));vmin=min(p_roi(2,:));vmax=max(p_roi(2,:));
options.KLims=[umin,umax,vmin,vmax]; % bounds
C=TPSGenerateCenters(delta.nC,options.KLims+1e-3.*[-1,1,-1,1]);
delta.C=C;
delta.ir=1e-4; % internal smoothing 1e-4
delta.er=1e-4; % external smoothing 1e -4
% Precompute EpsilonLambda matrix
EpsilonLambda=TPSEpsilonLambda(delta.C,delta.ir);
delta.EpsilonLambda=EpsilonLambda;
% Get delta parameters from correspondences
L=TPSWfromfeatures(p_roi,P_roi,delta.C,delta.er,delta.ir,delta.EpsilonLambda);
delta.L=L;
 [J_delta,~]=TPSWarpDiff(p_roi,delta.L,delta.C,delta.ir,delta.EpsilonLambda);

options.eta.er=e; % eta external smoothing %500 for book
options.eta.nC=5; % number of control centers
options.phi.er=p; % phi external smoothing % 0.0001 for book dnt remember
options.phi.nC=5; % number of control centers
options.maxiter=50; % max number of iterations with method='ReIso'
options.isoer=1e8; % lagrangian parameter for nonlinear refinement method='ReIso'
options.method='AnIso'; % 'ReIso' for nonlinear refinement
options.delta=delta;
 %options.verbose=1;
out=SfT(p_roi,q,options); % Shape from Template main function
phi=out.phi;
% Get 3D coordinates of p
[J_phi,Qw]=TPSWarpDiff(p_roi,phi.L,phi.C,phi.ir,phi.EpsilonLambda);
Q=Qw*scale;
