function [K] = getTangentOperator(geometry,material,D,strain,type)
array2tensor=@(strain) [strain(1) strain(3); strain(3) strain(2)];
switch type
    case 'fullY'
TangentOperator=@(G,K,D,trStrain) 2*G*(1-D)*[2/3 -1/3 0; -1/3 2/3 0; 0 0 1/2]...
    +(1-heaviside(trStrain))*K*[1 1 0; 1 1 0;0 0 0]...
    +heaviside(trStrain)*(1-D)*K*[1 1 0; 1 1 0;0 0 0];
    case 'H'
TangentOperator=@(G,K,D,trStrain) 2*G*[2/3 -1/3 0; -1/3 2/3 0; 0 0 1/2]...
    +(1-heaviside(trStrain))*K*[1 1 0; 1 1 0;0 0 0]...
    +heaviside(trStrain)*(1-D)*K*[1 1 0; 1 1 0;0 0 0];
end
K=zeros(geometry.maxEdof,geometry.maxEdof);
for j=1:geometry.Nele
%%approach 1
    Ke=reshape(permute(geometry.B(:,:,:,j),[2,1,3]),8,12)*...
    blkdiag(...
    TangentOperator(material.G,material.K,D(geometry.EGPconn(j,1)),trace(array2tensor(strain(geometry.EGPconn(j,1),:))))*geometry.wp(1),... %operator @ Gp1
    TangentOperator(material.G,material.K,D(geometry.EGPconn(j,2)),trace(array2tensor(strain(geometry.EGPconn(j,2),:))))*geometry.wp(2),... %operator @ Gp2
    TangentOperator(material.G,material.K,D(geometry.EGPconn(j,3)),trace(array2tensor(strain(geometry.EGPconn(j,3),:))))*geometry.wp(3),... %operator @ Gp3
    TangentOperator(material.G,material.K,D(geometry.EGPconn(j,4)),trace(array2tensor(strain(geometry.EGPconn(j,4),:))))*geometry.wp(4)...  %operator @ Gp4
    )*reshape(permute(geometry.B(:,:,:,j),[2,1,3]),8,12)'*geometry.detJ(j)*geometry.t;

%%apporach 2
% Ke=zeros(8,8);
% for i=1:geometry.NGPperElement %1:4    
%     Ke=Ke+transpose(geometry.B(:,:,i,j))*TangentOperator(material.G,material.K,D(geometry.EGPconn(j,i)),trace(array2tensor(strain(geometry.EGPconn(j,i),:))))*geometry.B(:,:,i,j)*geometry.wp(i);
% end
% Ke=Ke*geometry.detJ(j)*geometry.t;


%assemble global K matrix 
%K=assem(geometry.Edof(j,:),zeros(geometry.maxEdof,geometry.maxEdof),Ke)+K;
K=assem(geometry.Edof(j,:),K,Ke);
end
end

