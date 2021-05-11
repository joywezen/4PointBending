function [geometry]= setObject(ObjectID,geometry)

switch ObjectID
    case 1
   geometry.Coord=[0,0;
                1,0;
                2,0;
                3,0;
                0,1;
                1,1;
                2,1;
                3,1;
]   ;
geometry.Dof=[ 1  2;
      3  4;
      5  6;
      7  8;
      9  10;
      11 12
      13 14
      15 16];
geometry.Econn=[ 1 2 6 5;
        2 3 7 6;
        3 4 8 7;];

  %boundary conditions
    geometry.fixednodes=[1,5];
    geometry.loadnodes=[4,8];
    geometry.directions=geometry.loadDirection;
    [geometry.Xgrid,geometry.Ygrid]=meshgrid(linspace(0,3,4),linspace(0,1,2));
    case 2 %left fixed shell
        findGp= @(element) sort(reshape(repmat((element-1).*4,1,4)+repmat((1:4),size(element,1),1),size(element,1)*4,1));
           [geometry]=genMeshBar(geometry);
             
             geometry.fixednodes=findNode(geometry,geometry.fixCoord);%[(size(geometry.nodeSeq,2)-1)/5,size(geometry.nodeSeq,2)-(size(geometry.nodeSeq,2)-1)/5+1];
             geometry.loadnodes=findNode(geometry,geometry.loadCoord);
          
             
             %geometry.fixedEle=sort([geometry.eleSeq(geometry.NeleH,geometry.NeleL/5+1:geometry.NeleL/5*4)';reshape((geometry.eleSeq(1:geometry.NeleH,1:geometry.NeleL/5)),geometry.NeleH*geometry.NeleL/5,1);...
             %                  reshape((geometry.eleSeq(1:geometry.NeleH,geometry.NeleL/5*4+1:geometry.NeleL)),geometry.NeleH*geometry.NeleL/5,1)]);
                
    
    case 3
    [~,~,nodes,coordinates]=callCircularPlate();
    geometry.Coord=coordinates;  
    geometry.Dof=[1:2:numel(coordinates)-1;2:2:numel(coordinates)]';
    geometry.Econn=nodes;
    
    %boundary conditions
    geometry.fixednodes=[1,2,12,22, 32, 42 ,52,62];
    geometry.loadnodes=[11:10:61];
    geometry.directions=[1];
    
    %geometry.bc.dof=reshape(geometry.Dof([1,2,12,22, 32, 42 ,52,62],:),1,numel(geometry.Dof([1,2,12,22, 32, 42 ,52,62],:)));
%geometry.bc.U=zeros(1,numel(geometry.Dof([1,2,12,22, 32, 42 ,52,62],:)));

% geometry.Edof=[transpose(1:size(geometry.Econn,1)),transpose(reshape(geometry.Dof(geometry.Econn',:)',8,size(geometry.Econn,1)))];
%     geometry.maxEdof=max(max(geometry.Edof(2:end,:)));%maximum 
   
end
 geometry.Edof=[transpose(1:size(geometry.Econn,1)),transpose(reshape(geometry.Dof(geometry.Econn',:)',8,size(geometry.Econn,1)))];
    geometry.maxEdof=max(max(geometry.Edof(2:end,:)));%maximum 
    
geometry.bc.dof=reshape(geometry.Dof(geometry.fixednodes,:),1,numel(geometry.Dof(geometry.fixednodes,:)));  
geometry.bc.dof(2:end/2)=[];
geometry.bc.U=zeros(1,numel(geometry.Dof(geometry.fixednodes,:)));
geometry.bc.U(2:end/2)=[];
geometry.freedof=(1:geometry.maxEdof);
geometry.freedof(geometry.bc.dof)=[];

geometry.ir=2;
geometry.NGPperElement=geometry.ir^geometry.ProblemDimension;
geometry.NGPtotal=geometry.NGPperElement*geometry.Nele;
%  for i=1:size(geometry.Econn,1)
%  geometry.CP(:,i)=mean(geometry.Coord(geometry.Econn(i,:)',:));
%  end
 geometry.CoordEle=(geometry.Coord(geometry.Econn',:));

 %##method 1 to get Edof
 %geometry.Edof=[1 reshape(geometry.Dof(geometry.Econn(1,:),:)',1,8)
  %     2 reshape(geometry.Dof(geometry.Econn(2,:),:)',1,8)];
 %##method 2 to get Edof
 
  % ----- Element properties, topology and coordinates -----
geometry.nen=4;
[geometry.Ex,geometry.Ey]=coordxtr(geometry.Edof,geometry.Coord,geometry.Dof,geometry.nen);
%clf; 
%eldraw2(geometry.Ex,geometry.Ey,[1 2 1],geometry.Edof(:,1));
geometry.t=0.1;%thickness
             [geometry]=getGp(geometry);
           %  geometry.NoDamageGp=findGp(geometry.fixedEle);
             geometry.DamageGp=transpose(1:geometry.NGPtotal); %%%%%%%%%
             %geometry.DamageGp(geometry.NoDamageGp)=[];
             geometry.DamageGp=transpose(1:geometry.NGPtotal);
             geometry.EGPconn=transpose([((1:geometry.Nele)-1)*4+1;((1:geometry.Nele)-1)*4+2;...
                                         ((1:geometry.Nele)-1)*4+3;((1:geometry.Nele)-1)*4+4]);
end