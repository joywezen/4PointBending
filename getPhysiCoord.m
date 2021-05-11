function [PhysiCoordx,PhysiCoordy]=getPhysiCoord(geometry,localx,localy)
  if nargin==1 %get GaussPoint physical coordinate
  xsi=geometry.GP(:,1);
  eta=geometry.GP(:,2);
  else % get other physical coordinate 
  xsi=localx;
  eta=localy;
  end
  N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
  N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;
  
  %N=[ N1 N2 N3 N4  for GP1
  %    N1 N2 N3 N4  for GP2...
  %...
  %    N1 N2 N3 N4] for GP4
  
  PhysiCoordx=N*geometry.Ex';% n GPs per element * n of elements
  PhysiCoordy=N*geometry.Ey';% n GPs per element * n of elements
end