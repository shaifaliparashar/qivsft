function s=find_scale(p_roi,P_roi,Q_roi)
Mdl = createns(p_roi,'NSMethod','kdtree','Distance','euclidean');
%Mdl = createns(P_roi,'NSMethod','kdtree','Distance','euclidean');
% scale = A (p dist)\B (q_dist) 
A=[];
B=[];
for i=1:length(p_roi)
[Idx,D] = knnsearch(Mdl,p_roi(i,:),'Distance','euclidean','k',21);
p=repmat(P_roi(i,:),20,1);
q=repmat(Q_roi(i,:),20,1);
A= [A; sqrt(sum((p-P_roi(Idx(2:21),:)).^2,2))];
B= [B; sqrt(sum((q-Q_roi(Idx(2:21),:)).^2,2))];
end

s= A\B;