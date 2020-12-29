function [T,P_dis,Pest,L,idd,idd2,lambda,C,EpsilonLambda,T_G]= Init_volume(P_undis,P_dis,P_roi,Q_roi,alpha,d,th,muu1,lambda,Q_G,J_delta,J_phi,init)
                                                          %Init_volume(template',gth_liver',temp_3d_pts,ip_3d_est,alpha,d,th,1e-4,1e-4,ip_3d_pts,J_delta,J_phi,init);
%STEP 1: find points on surface
P_undis=P_undis*alpha;% make T on the same scale as the ground truth
P_roi=P_roi*alpha;
T=[P_undis;P_roi']; % template/
T=double(unique(T,'rows'));
DT = delaunayTriangulation(T(:,1),T(:,2),T(:,3));
T_G=[P_dis;Q_G'];

x= min(T(:,1))-1:d:max(T(:,1))+1;
y= min(T(:,2))-1:d:max(T(:,2))+1;
z= min(T(:,3))-1:d:max(T(:,3))+1;
yy=repmat(y',1,length(x));
yy=yy';
yy=yy(:);
zz=repmat(z',1,length(yy));
zz=zz';
zz=zz(:);
Pn= [repmat(x',length(y)*length(z),1),repmat(yy,length(z),1),zz];
ti= pointLocation(DT,Pn);
id=1:length(ti);
%Remove Pn outside the mesh
id(isnan(ti)==1)=[];
Padded=Pn(id,:);
T=[T;Padded];
T=double(unique(T,'rows'));
DT = delaunayTriangulation(T(:,1),T(:,2),T(:,3));
[~,idd]=ismember(P_undis,T,'rows');
[~,idd2]=ismember(P_roi',T,'rows');


%STEP 2: Propagate into the volume
DTS = struct('Points', DT.Points, 'ConnectivityList', DT.ConnectivityList);
for i=1:length(DT.ConnectivityList)
    t=DT.ConnectivityList(i,:);
    a=T(t(1),:);
    b=T(t(2),:);
    c=T(t(3),:);
    d=T(t(4),:);
    dis(1)=sqrt(sum((T(t(1),:)-T(t(2),:)).^2));
    dis(2)=sqrt(sum((T(t(1),:)-T(t(3),:)).^2));
    dis(3)=sqrt(sum((T(t(1),:)-T(t(4),:)).^2));
    dis(4)=sqrt(sum((T(t(2),:)-T(t(3),:)).^2));
    dis(5)=sqrt(sum((T(t(2),:)-T(t(4),:)).^2));
    dis(6)=sqrt(sum((T(t(3),:)-T(t(4),:)).^2));
    if (max(dis)-min(dis))> 1
        DTS.ConnectivityList(i,:)=[0 0 0 0];
    end
end

is = 1:length(DT.ConnectivityList);
is(sum(DTS.ConnectivityList,2)==0)=[];
DTS = struct('Points', DTS.Points, 'ConnectivityList', DTS.ConnectivityList(is,:));
cnt =1;
T_new=zeros(size(T));
[~ ,v1]=ismember(T,P_roi','rows');
idx=1:length(T);
idx(v1==0)=[];
if init <2
    T_new(idx,:)= Q_roi(:,v1(idx))';
end
if init == 0
    % th=.05;
    while cnt >0
        cnt =0;
        tr=[]; %list of traversed tetrahedrons
        idx=[];% list of tetrahedrons to be updated
        %     find tetrahedrons whose all pts are known
        V=DT.ConnectivityList;
        for i=1:length(DT.ConnectivityList)
            for j=1:4
                V(i,j)= sum(abs(T_new(DT.ConnectivityList(i,j),:)))>0;
            end
            if sum(V(i,:))>3
                tr=[tr,i];
            end
            if sum(V(i,:))==3
                idx=[idx,i];
            end
        end
        
        if ~isempty(idx)
            for i=1: length(idx)
                id =idx(i);
                [~,o]=sort(V(id,:));
                %             pts in template
                a1=T(DT.ConnectivityList(id,o(2)),:);
                b1=T(DT.ConnectivityList(id,o(3)),:);
                c1=T(DT.ConnectivityList(id,o(4)),:);
                %             pts in deformation
                a=T_new(DT.ConnectivityList(id,o(2)),:);
                b=T_new(DT.ConnectivityList(id,o(3)),:);
                c=T_new(DT.ConnectivityList(id,o(4)),:);
                d1=T(DT.ConnectivityList(id,o(1)),:);
                d(1)=norm(a1-d1);
                d(2)=norm(b1-d1);
                d(3)=norm(c1-d1);
                R=[a1;b1;c1]\[a;b;c];
                dnew=d1*R;
                dn(1)=norm(a-dnew);
                dn(2)=norm(b-dnew);
                dn(3)=norm(c-dnew);
                if abs(mean(dn)-mean(d)) < th && v1(DT.ConnectivityList(id,o(1)))==0
                    
                    T_new(DT.ConnectivityList(id,o(1)),:)=dnew;
                    cnt=cnt+1;
                else
                    idx(i)=0;
                    
                end
            end
            
        end
    end
    
    
    % plot3(T(:,1),T(:,2),T(:,3),'*b');
    s=sum(abs(T_new),2);
    s_i=1:length(T_new);
    s_i(s==0)=[];
    
    n=0;
    for i=1:length(T_new)
        if sum(T_new(i,:)==0)
            T_new(i,:)=T(i,:);
            n=n+1;
        end
    end
    %    n
    [~,id2]=ismember(Q_roi',T_new,'rows');
end

T=T';
T_new=T_new';


KLims= [min(T(1,:)) max(T(1,:)) min(T(2,:)) max(T(2,:)) min(T(3,:)) max(T(3,:)) ];
C=TPSGenerateCenters(3,KLims);
%
EpsilonLambda=TPSEpsilonLambda(C,lambda);
% prev s_i
if init == 0 %Greedy Initialisation
    [L,ZTZ]=TPSWfromfeatures(T(:,s_i),T_new(:,s_i),C,muu1,lambda,EpsilonLambda);
    
    [~,Pest]=TPSWarpDiff(T,L,C,lambda,EpsilonLambda);
    save('Pest.mat','Pest');
    save('idd.mat','idd');
    save('id2.mat','id2');
    save('T.mat','T');
    save('C.mat','C');
    save('EpsilonLambda.mat','EpsilonLambda');
    save('ZTZ.mat','ZTZ');
    
else if init ==1 % SMOOTH INITIALISATION
        [~,id2]=ismember(Q_roi',T_new','rows');
        [L,ZTZ]=TPSWfromfeatures([T(:,id2)],[T_new(:,id2)],C,muu1,lambda,EpsilonLambda);% initialise L
        [~,Pest]=TPSWarpDiff(T,L,C,lambda,EpsilonLambda);
        
        save('Pest.mat','Pest');
        save('idd.mat','idd');
        save('id2.mat','id2');
        
        save('T.mat','T');
        save('C.mat','C');
        save('EpsilonLambda.mat','EpsilonLambda');
        save('ZTZ.mat','ZTZ');
    end
end


