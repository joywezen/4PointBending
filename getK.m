function [K]=getK(geometry,material,D)
Ex=geometry.Ex;
%Ey=geometry.Ey;
t=geometry.t;
B=geometry.B;
%C=material.C;
Edof=geometry.Edof;
%maxEdof=geometry.maxEdof;
K=zeros(geometry.maxEdof,geometry.maxEdof);
%Ke=zeros(8,8,size(Ex,1));
% materialE=material.E;
% materialnv=material.nv;
% geometrywp1=geometry.wp(1);
% geometrywp2=geometry.wp(2);
% geometrywp3=geometry.wp(3);
% geometrywp4=geometry.wp(4);
% geometrydetJ=geometry.detJ;
for j=1:size(Ex,1)
%[Ke,~,B(:,:,:,j)]=planre(Ex(j,:),Ey(j,:),Ep,C,[],x((j-1)*4+(1:4))); %  damage varjables for each GP
%a=(Ex(j,3)-Ex(j,1))/2;  
%b=(Ey(j,3)-Ey(j,1))/2;

% for i=1:size(Ex,2)
% Ke=Ke+B(:,:,i,j)'*C*(1-D(i,j))*B(:,:,i,j)*a*b*t;%%see line 8 of inFEDma
% end
%tic
Ke=reshape(permute(B(:,:,:,j),[2,1,3]),8,12)*...
    blkdiag(...
    hooke(2,material.E,material.nv)*(1-D(geometry.EGPconn(j,1)))*geometry.wp(1),...
    hooke(2,material.E,material.nv)*(1-D(geometry.EGPconn(j,2)))*geometry.wp(2),...
    hooke(2,material.E,material.nv)*(1-D(geometry.EGPconn(j,3)))*geometry.wp(3),...
    hooke(2,material.E,material.nv)*(1-D(geometry.EGPconn(j,4)))*geometry.wp(4))*...
    reshape(permute(B(:,:,:,j),[2,1,3]),8,12)'*geometry.detJ(j)*t;
%toc
% Ke=reshape(permute(B(:,:,:,j),[2,1,3]),8,12)*...
%     blkdiag(hooke(2,materialE,materialnv)*(1-D((j-1)*4+1))*geometrywp1,hooke(2,materialE,materialnv)*(1-D((j-1)*4+2))*geometrywp2,...
%     hooke(2,materialE,materialnv)*(1-D((j-1)*4+3))*geometrywp3,hooke(2,materialE,materialnv)*(1-D((j-1)*4+4))*geometrywp4)*...
%     reshape(permute(B(:,:,:,j),[2,1,3]),8,12)'*geometrydetJ(j)*t;
%tic
%K=assem(Edof(j,:),zeros(maxEdof,maxEdof),Ke)+K;
K=assem(Edof(j,:),K,Ke);
%toc
end
%K=sum(K,3);
end