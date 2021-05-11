function [geometry]=getBmat(geometry)

B=zeros(3,8,size(geometry.Econn,2),size(geometry.Econn,1));
detJ=zeros(1,size(geometry.Econn,1));
ex=geometry.Ex;
ey=geometry.Ey;
% [xgp,ygp]=getGp;
% a=zeros(size(ey,1),1);
% b=zeros(size(ey,1),1);
%      for j=1:size(ey,1)
%      
%      a(j)=(ex(j,3)-ex(j,1))/2;  b(j)=(ey(j,3)-ey(j,1))/2;
%      
%      for i=1:size(ex,2)   
%        x=xgp(i)*a(j); y=ygp(i)*b(j);
%     
%        B(:,:,i,j)=[-(b(j)-y)    0     b(j)-y     0   b(j)+y   0  -(b(j)+y)   0  ;
%              0   -(a(j)-x)    0   -(a(j)+x)  0   a(j)+x    0    a(j)-x ;
%           -(a(j)-x) -(b(j)-y) -(a(j)+x)   b(j)-y  a(j)+x  b(j)+y   a(j)-x -(b(j)+y)]/(4*a(j)*b(j));
%      end
%      
%      end
  ngp=geometry.ir*geometry.ir;
  %w=geometry.w;
  gp=geometry.GP;
  
  %wp=w(:,1).*w(:,2);
  xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;
  
  %--------- shape functions -----------------------------------
  N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
  N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;

  dNr(1:2:r2,1)=-(1-eta)/4;     dNr(1:2:r2,2)= (1-eta)/4;
  dNr(1:2:r2,3)= (1+eta)/4;     dNr(1:2:r2,4)=-(1+eta)/4;
  dNr(2:2:r2+1,1)=-(1-xsi)/4;   dNr(2:2:r2+1,2)=-(1+xsi)/4;
  dNr(2:2:r2+1,3)= (1+xsi)/4;   dNr(2:2:r2+1,4)= (1-xsi)/4;
  
  dN2dxsi=[-(1)/4 (1)/4 (1)/4 -(1)/4];
  dN2deta=[-(1)/4 -(1)/4 (1)/4 (1)/4];
  for j=1:size(ey,1)
  JT=[dN2dxsi;dN2deta] * [ex(j,:);ey(j,:)]';
  detJ(j)=det(JT);
  %--------- plane stress --------------------------------------
  for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      
      %JTinv=inv(JT(indx,:));
      %dNx=JTinv*dNr(indx,:);
      
      dNx=JT\dNr(indx,:);

      B(1,1:2:8-1,i,j)=dNx(1,:);
      B(2,2:2:8,i,j)  =dNx(2,:);
      B(3,1:2:8-1,i,j)=dNx(2,:);
      B(3,2:2:8,i,j)  =dNx(1,:);

      %N2(1,1:2:8-1)=N(i,:);
      %N2(2,2:2:8)  =N(i,:);

     % Ke=Ke+B'*Dm*B*detJ*wp(i)*t;
     % fe=fe+N2'*b*detJ*wp(i)*t;
  end
  
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
  end
  geometry.B=B;
  geometry.detJ=detJ;
  
end