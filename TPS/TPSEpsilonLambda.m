function Epsilon_lambda=TPSEpsilonLambda(C,lambda)
l=size(C,1);
Ct=[C,ones(l,1)];
Klambda=lambda.*eye(l);
%obtaining K_lambda
for i=[1:l]
    for j=[1:l]
        if(i~=j)
            d=sum((C(i,:)-C(j,:)).^2);
            Klambda(i,j)=TPSrho(d);
        end
    end
    end
%A=inv(Ct'*inv(Klambda)*Ct)*Ct'*inv(Klambda);
A=((Ct'*(Klambda\Ct))\Ct')/Klambda;
Epsilon_lambda=[(Klambda)\(eye(l)-Ct*A);A];