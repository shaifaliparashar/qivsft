function Pest= Refinement(P_dis,T,L,muu,rho,lambda,C,EpsilonLambda, ref)
x=L(:);
save('muu.mat','muu');
save('rho.mat','rho');
opt = optimset('Jacobian','on','Display','iter','MaxIter',200, 'DerivativeCheck','off');
if ref== 2 % L2 norm
[Lnew,J]=lsqnonlin(@myfunED3syn,x,[],[],opt);
else if ref==1 % L1 norm
        [Lnew,J]=lsqnonlin(@myfunED3syn2,x,[],[],opt);
    end
end 
L=Lnew;
L=reshape(L,size(C,1),3);
[~,Pest]=TPSWarpDiff(T,L,C,lambda,EpsilonLambda);
% % 
% figure;hold on
% plot3(Pest(1,:),Pest(2,:),Pest(3,:),'o','LineWidth',2,'Color',[0.6,0,0]);
% plot3(P_dis(:,1),P_dis(:,2),P_dis(:,3),'o','LineWidth',2,'Color',[0,0.6,0]);
% legend('Estimated Points','Ground Truth')
% axis equal
% hold off
