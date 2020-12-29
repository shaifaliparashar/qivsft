function [P_roi,Q_roi,p_roi,q_roi,q]=find3d(p_roi,q_roi,Undis,Dis,q)

P_dis=Dis.P.vertexPos; % 3D points on the deformed model

P_undis= Undis.P.vertexPos; % 3D points on template
Q_roi=[]; % Ground truth
P_roi=[]; % corresponding 3D points
idx=[];

for i=1:length(p_roi)
  
    bp=Undis.P.B(round(p_roi(2,i)),round(p_roi(1,i)),:);% barycentrics of p_roi
    bpQ=Dis.P.B(round(q_roi(2,i)),round(q_roi(1,i)),:);% barycentrics of q_roi
    % bp=[facen,alpha1,alpha2,alpha3]
    % get triangle points'
    if bp(1)<=0 || bpQ(1)<=0
        idx=[idx;i];
         
    else
        % for P_roi
        facen=Undis.P.faces(bp(1),:);
        P1=P_undis(facen(1),:)';
        P2=P_undis(facen(2),:)';
        P3=P_undis(facen(3),:)';
        P_roi=[P_roi,bp(2).*P1+bp(3).*P2+bp(4).*P3];
         % for Q_roi
        facenQ=Dis.P.faces(bpQ(1),:);
        P4=P_dis(facenQ(1),:)';
        P5=P_dis(facenQ(2),:)';
        P6=P_dis(facenQ(3),:)';
        Q_roi=[Q_roi,bpQ(2).*P4+bpQ(3).*P5+bpQ(4).*P6];
        
    end
end
p_roi(:,idx)=[];
q_roi(:,idx)=[];
q(:,idx)=[];

