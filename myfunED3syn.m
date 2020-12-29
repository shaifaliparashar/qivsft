function [F,J] = myfunED3syn(x)

load('C.mat')
load('T.mat')
load('wP.mat')
load('ZTZ.mat')
load('EpsilonLambda.mat')
load('id2.mat')
load('rho.mat')
load('muu.mat')


lambda=1e-4;
[U, S, V]= svd(ZTZ);
D=sqrt(S);
x=reshape(x,size(U,1),3);
l=size(wP,2);
%keyboard
%[dP, Pest]=TPSWarpDiff(T,x,C,lambda,EpsilonLambda);
[M1, M2]=TPSCoeffDiff(T,x,C,lambda,EpsilonLambda);
Pest= M1*x;
Pest=Pest';
dP=M2*x;

Fd=(wP-[Pest(1,id2)./Pest(3,id2);Pest(2,id2)./Pest(3,id2)])/sqrt(l);

%jacobian for data term
Jd=[];
for i=1:l
    ti=[-1/Pest(3,id2(i)) 0 Pest(1,id2(i))/(Pest(3,id2(i)).^2); 0 -1/Pest(3,id2(i)) Pest(2,id2(i))/(Pest(3,id2(i)).^2)]*blkdiag(M1(id2(i),:),M1(id2(i),:),M1(id2(i),:));
    Jd=[Jd;ti]; 
end
% figure; hold on
% plot(wP(1,:),wP(2,:),'*r');
% plot(Pest(1,id2)./Pest(3,id2),Pest(2,id2)./Pest(3,id2),'*g');
% hold off
% axis equal
Jd=Jd/sqrt(l);
Fs= D*V'*x;
Js=D*V';

Js=blkdiag(Js,Js,Js);

Ji=[];
Fi=[];

% Pi=Pest(:,idd);
% idx=[idd, idd+1,idd+2];
% idx=idx';
% idx=idx(:);
% dP=dP(idx,:);
% M2=M2(idx,:);

for i=1:size(Pest,2)
    dPe=dP((i-1)*3 +1:(i-1)*3+ 3,:)';
    dPe= dPe'*dPe - eye(3);
    dPe=[dPe(1,1);dPe(1,2);dPe(1,3);dPe(2,2);dPe(2,3);dPe(3,3)]; % upper
    %     
    % dPe(1,1)=(M2(1,:)*x)*(M2(1,:)*x)'
    M2i=M2((i-1)*3 + 1:(i)*3,:);
    % dPde(1,1)
    Jii=[2*(M2i(1,:)*x(:,1))*M2i(1,:),2*(M2i(1,:)*x(:,2))*M2i(1,:),2*(M2i(1,:)*x(:,3))*M2i(1,:)];
    % dPde(1,2)
    Jii=[Jii;[(M2i(1,:)*x(:,1))*M2i(2,:)+(M2i(2,:)*x(:,1))*M2i(1,:),(M2i(1,:)*x(:,2))*M2i(2,:)+(M2i(2,:)*x(:,2))*M2i(1,:),(M2i(1,:)*x(:,3))*M2i(2,:)+(M2i(2,:)*x(:,3))*M2i(1,:)]];
     % dPde(1,3)
    Jii=[Jii;[(M2i(1,:)*x(:,1))*M2i(3,:)+(M2i(3,:)*x(:,1))*M2i(1,:),(M2i(1,:)*x(:,2))*M2i(3,:)+(M2i(3,:)*x(:,2))*M2i(1,:),(M2i(1,:)*x(:,3))*M2i(3,:)+(M2i(3,:)*x(:,3))*M2i(1,:)]];
     % dPde(2,2)
    Jii=[Jii;[(M2i(2,:)*x(:,1))*M2i(2,:)+(M2i(2,:)*x(:,1))*M2i(2,:),(M2i(2,:)*x(:,2))*M2i(2,:)+(M2i(2,:)*x(:,2))*M2i(2,:),(M2i(2,:)*x(:,3))*M2i(2,:)+(M2i(2,:)*x(:,3))*M2i(2,:)]];    
     % dPde(2,3)
    Jii=[Jii;[(M2i(2,:)*x(:,1))*M2i(3,:)+(M2i(3,:)*x(:,1))*M2i(2,:),(M2i(2,:)*x(:,2))*M2i(3,:)+(M2i(3,:)*x(:,2))*M2i(2,:),(M2i(2,:)*x(:,3))*M2i(3,:)+(M2i(3,:)*x(:,3))*M2i(2,:)]];    
     % dPde(3,3)
    Jii=[Jii;[(M2i(3,:)*x(:,1))*M2i(3,:)+(M2i(3,:)*x(:,1))*M2i(3,:),(M2i(3,:)*x(:,2))*M2i(3,:)+(M2i(3,:)*x(:,2))*M2i(3,:),(M2i(3,:)*x(:,3))*M2i(3,:)+(M2i(3,:)*x(:,3))*M2i(3,:)]];    
    Ji=[Ji;Jii];
    %triangular
%     %dPe=[dPe(1,1);dPe(1,2);dPe(1,3);dPe(2,2);dPe(2,3);dPe(3,3)]; %lower triangular
%     t(1:3,:)=(2*M2((i-1)*3 + 1,:)'*M2((i-1)*3 + 1,:)*x)';
%     %t(1:3,:)=x'*M2((i-1)*3 + 1,:)'*M2((i-1)*3 + 1,:) + (M2((i-1)*3 + 1,:)'*M2((i-1)*3 + 1,:)*x)';
%     Ji((i-1)*6 + 1,:)=[t(1,:),t(2,:),t(3,:)];
%     t(4:6,:)=x'*M2((i-1)*3 + 1,:)'*M2((i-1)*3 + 2,:) + (M2((i-1)*3 + 1,:)'*M2((i-1)*3 + 2,:)*x)';
%     Ji((i-1)*6 + 2,:)=[t(4,:),t(5,:),t(6,:)];
%     t(7:9,:)=(2*M2((i-1)*3 + 2,:)'*M2((i-1)*3 + 2,:)*x)';
%     %t(7:9,:)=x'*M2((i-1)*3 + 2,:)'*M2((i-1)*3 + 2,:) + (M2((i-1)*3 + 2,:)'*M2((i-1)*3 + 2,:)*x)';
%     Ji((i-1)*6 + 3,:)=[t(7,:),t(8,:),t(9,:)];
%     t(10:12,:)=x'*M2((i-1)*3 + 3,:)'*M2((i-1)*3 + 1,:) + (M2((i-1)*3 + 3,:)'*M2((i-1)*3 + 1,:)*x)';
%     Ji((i-1)*6 + 4,:)=[t(10,:),t(11,:),t(12,:)];
%     t(13:15,:)=x'*M2((i-1)*3 + 3,:)'*M2((i-1)*3 + 2,:) + (M2((i-1)*3 + 3,:)'*M2((i-1)*3 + 2,:)*x)';
%     Ji((i-1)*6 + 5,:)=[t(13,:),t(14,:),t(15,:)];
%     t(16:18,:)=(2*M2((i-1)*3 + 3,:)'*M2((i-1)*3 + 3,:)*x)';
%     %t(16:18,:)=x'*M2((i-1)*3 + 3,:)'*M2((i-1)*3 + 3,:) + (M2((i-1)*3 + 3,:)'*M2((i-1)*3 + 3,:)*x)';
%     Ji((i-1)*6 + 6,:)=[t(16,:),t(17,:),t(18,:)];
    Fi=[Fi;dPe];
end

%Fi=Fi';
% for i=1:size(P,2)
%     dPest= dP(:,i);
%     dPest=reshape(dPest,3,3);
%     dPest=dPest*dPest'-eye(3);
%     dPest=dPest(:);
%     fi=[fi,dPest(1),dPest(4),dPest(5),dPest(7),dPest(8),dPest(9)];
% end
 %a=sqrt(Fd(:)'*Fd(:))
% b=sqrt(Fs(:)'*Fs(:))
% c=sqrt(Fi*Fi')
F=[ Fd(:);muu*Fs(:);rho*Fi];
J=[Jd;muu*Js;rho*Ji];
save('x.mat','x');
end