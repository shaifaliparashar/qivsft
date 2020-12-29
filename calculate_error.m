function [err,P_roi,Q_roi]= calculate_error(T,T_new,gp_roi,gq_roi,P_undis,Dis,GUndis,GDis,alpha)
%[err_ref,P_roi,Q_roi]=calculate_error(T',Pest',temp_back_pts,ip_back_pts,temp_ver,Input,Template_back,Input_back,alpha);
% gp_roi=gp_roi';
% gq_roi=gq_roi';
P_undis=P_undis*alpha;
P_dis=Dis.P.vertexPos; % 3D points on the deformed model

iid=[];
for i=1:length(P_undis)
    m=sum((repmat(P_undis(i,:),length(T),1)-T).^2,2);
    [~, ix]=min(m);
    iid=[iid,ix];
end
Q_roi=[]; % Ground truth
P_roi=[]; % corresponding 3D points


for i=1:length(gp_roi)
    
    bp=GUndis.P.B(round(gp_roi(i,2)),round(gp_roi(i,1)),:);% barycentrics of p_roi
    bpQ=GDis.P.B(round(gq_roi(i,2)),round(gq_roi(i,1)),:);% barycentrics of q_roi
    
    % get triangle points'
   
        % for P_roi
        if bp(1) >0 && bpQ(1)
        facen=GUndis.P.faces(bp(1),:);
        P1=T_new(iid(facen(1)),:)';
        P2=T_new(iid(facen(2)),:)';
        P3=T_new(iid(facen(3)),:)';
        P_roi=[P_roi,bp(2).*P1+bp(3).*P2+bp(4).*P3];
         % for Q_roi
        facenQ=GDis.P.faces(bpQ(1),:);
        P4=P_dis(facenQ(1),:)';
        P5=P_dis(facenQ(2),:)';
        P6=P_dis(facenQ(3),:)';
        Q_roi=[Q_roi,bpQ(2).*P4+bpQ(3).*P5+bpQ(4).*P6];
        end
    
end
err= sqrt(mean(sum((P_roi- Q_roi).^2)));% 

