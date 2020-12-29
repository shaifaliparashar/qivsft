% Example script for Volumetric Shape-from_Template
% By: Shaifali Parashar

addpath('./TPS');
addpath('./SfT');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% Input image and ground truth front
input=imread('input/IMG_1819.tif');
Input=load('input/IMG_1819_GT.mat','P');


% Template front
template=imread('input/IMG_1607.tif');
Template=load('input/IMG_1607_GT.mat','P');

% Template back
template_back=imread('input/IMG_1663.tif');
Template_back=load('input/IMG_1663_GT.mat','P');

% Input back
input_back=imread('input/IMG_1875.tif');
Input_back=load('input/IMG_1875_GT.mat','P');

% surface points front
temp = load('input/fixedPoints.mat','fixedPoints');
ip = load('input/movingPoints.mat','movingPoints');
temp_surf_pts = temp.fixedPoints;
ip_surf_pts = ip.movingPoints;

% surface points back
temp = load('input/fixedPoints1.mat','fixedPoints');
ip = load('input/movingPoints1.mat','movingPoints');
temp_back_pts = temp.fixedPoints;
ip_back_pts = ip.movingPoints;

%%%%%%%%%%%% read parameters %%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('input/input_greedy_L1.txt');
a=textscan(fid,'%s',11);
fclose(fid);

b= a{1};
len= str2num(b{1});     % length of the object
scale= str2num(b{2});   % scale (NOT REQUIRED)
eta = str2num(b{3});    % eta SFT should be tuned
phi = str2num(b{4});    % phi SFT should be tuned
d = str2num(b{5});      % for sampling while initialising volume (0.1 CORRESPONDS TO SHOULD BE 3-4 MM)
th = str2num(b{6});     % threshold for greedy init
muu = str2num(b{7});    % refinement % SMOOTHNESS SHOULD BE SMALL
rho = str2num(b{8});    % refinement % ISOMETRY SHOULD BE SMALL
al = str2num(b{9});     % NOT REQUIRED
init= str2num(b{10});   % 0 for greedy (GREEDY FOR ACCURACY BUT SLOW)
ref=str2num(b{11});     % 1 for L1, 2 for L2


%% %%%%%%%%%%%%%%%%%%%%%%% Vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_ver = Template.P.vertexPos;
in_ver = Input.P.vertexPos;
temp_back_ver= Template_back.P.vertexPos;
in_back_ver= Input_back.P.vertexPos;
load Calib_Results KK

%%%%%%%%%%%%%% Points projected by the camera %%%%%%%%%%%%%%%%%%%%%%%

size(ip_surf_pts)

ip_2d_est = KK\[ip_surf_pts';ones(1,length(ip_surf_pts))];
ip_2d_est(3,:)=[];

%% %%%%%%%%%%% render 2D image points to 3D %%%%%%%%%%%%%%%%%%%%%%%%% (COMMENTS)
% find 3D of the surface points SFT
[temp_3d_pts,ip_3d_pts,temp_surf_pts,ip_surf_pts,ip_2d_est]=find3d(temp_surf_pts',ip_surf_pts',Template,Input,ip_2d_est);                                                                  
temp_3d_pts= double(temp_3d_pts); % TEMPLTE 3d
ip_3d_pts=double(ip_3d_pts); % GROUND TRUTH

%% %%%%%%%%%%%%%%%%%%%%%% find scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=find_scale(temp_surf_pts',temp_3d_pts',ip_3d_pts');


%%%%%%%%%%%%%%%%%%%% Surface Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ip_3d_est,J_delta,J_phi]=create_template(temp_surf_pts, temp_3d_pts, ip_2d_est, ip_3d_pts, eta,phi,alpha,KK);
scale = max(max(Input.P.vertexPos) - min(Input.P.vertexPos));
sft_err = sqrt(mean(sum((ip_3d_pts - ip_3d_est).^2)))*10*len/scale;
    
figure; title('Surface Reconstruction'); hold on;
plot3(ip_3d_est(1,:),ip_3d_est(2,:),ip_3d_est(3,:),'o','LineWidth',2,'Color',[0.6,0,0]);
plot3(ip_3d_pts(1,:),ip_3d_pts(2,:),ip_3d_pts(3,:),'o','LineWidth',2,'Color',[0,0.6,0]);
legend('Reconstruction','Ground Truth'); axis equal; hold off

%%%%%%%%%%%%%%%%%%%% Volume Propagation %%%%%%%%%%%%%%%%%%%%%%%%%
wP=ip_2d_est;
save('wP.mat','wP');
% tube distorted 1 and 2
[T,P_dis,T_new,L,idd,idd2,lambda,C,EpsilonLambda]= Init_volume(temp_ver,in_ver,temp_3d_pts,ip_3d_est,alpha,d,th,1e-4,1e-4,ip_3d_pts,J_delta,J_phi,init);

% show initialization
figure; title('Initialization'); hold on;
%plot3(T(1,:),T(2,:),T(3,:),'o','LineWidth',2,'Color',[0,0,0.6]);
plot3(T_new(1,:),T_new(2,:),T_new(3,:),'o','LineWidth',2,'Color',[0.6,0,0]);
plot3(P_dis(:,1),P_dis(:,2),P_dis(:,3),'o','LineWidth',2,'Color',[0,0.6,0]);
legend('Reconstruction','Ground Truth'); axis equal; hold off;

% initialization back surface
[err_init, P_roi, Q_roi]=calculate_error(T',T_new',temp_back_pts,ip_back_pts,temp_ver,Input,Template_back,Input_back,alpha);
err_init= err_init*10*len/scale
figure; hold on; 
plot3(P_roi(1,:),P_roi(2,:),P_roi(3,:),'o','LineWidth',2,'Color',[0.6,0,0]);
plot3(Q_roi(1,:),Q_roi(2,:),Q_roi(3,:),'o','LineWidth',2,'Color',[0,0.6,0]);
axis equal; hold off;

% initializatio front surface
err_init_surf=10*len*sqrt(mean(sum((ip_3d_pts - T_new(:,idd2)).^2)))/scale
figure; hold on;
plot3(T_new(1,idd2),T_new(2,idd2),T_new(3,idd2),'o','LineWidth',2,'Color',[0.6,0,0]);
plot3(ip_3d_pts(1,:),ip_3d_pts(2,:),ip_3d_pts(3,:),'o','LineWidth',2,'Color',[0,0.6,0]);
legend('Reconstruction','Ground Truth'); axis equal; hold off;

%%%%%%%%%%%%%%%%%%%%%% Refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pest= Refinement(P_dis,T,L,muu,rho,lambda,C,EpsilonLambda,ref);

% [~,sca,~]=RegisterToGTH(Pest(:,idd2),ip_3d_pts)
% %Pest=Pest*sca;
figure; title('Refinement'); hold on; 
%plot3(T(1,:),T(2,:),T(3,:),'o','LineWidth',2,'Color',[0,0,0.6]);
plot3(Pest(1,:),Pest(2,:),Pest(3,:),'o','LineWidth',2,'Color',[0,0,0.6]);
plot3(P_dis(:,1),P_dis(:,2),P_dis(:,3),'o','LineWidth',2,'Color',[0,0.6,0]);
legend('Reconstruction','Ground Truth'); axis equal; hold off

[err_ref,P_roi,Q_roi]=calculate_error(T',Pest',temp_back_pts,ip_back_pts,temp_ver,Input,Template_back,Input_back,alpha);
err_ref= err_ref*10*len/scale

%refinement front
figure; hold on; 
plot3(P_roi(1,:),P_roi(2,:),P_roi(3,:),'o','LineWidth',2,'Color',[0,0,0.6]);
plot3(Q_roi(1,:),Q_roi(2,:),Q_roi(3,:),'o','LineWidth',2,'Color',[0,0.6,0]);
axis equal; hold off;
err_ref_surf = 10*len*sqrt(mean(sum((ip_3d_pts - Pest(:,idd2)).^2)))/scale

%refinement back
figure; hold on
plot3(Pest(1,idd2),Pest(2,idd2),Pest(3,idd2),'o','LineWidth',2,'Color',[0,0,0.6]);
plot3(ip_3d_pts(1,:),ip_3d_pts(2,:),ip_3d_pts(3,:),'o','LineWidth',2,'Color',[0,0.6,0]);
%legend('Reconstruction','Ground Truth')
axis equal; hold off




