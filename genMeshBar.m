function [geometry]=genMeshBar(geometry)
L=geometry.L;
H=geometry.H;
L0=geometry.L0;
H0=geometry.H0;
NeleL=geometry.NeleL;
NeleH=geometry.NeleH;
[Xgrid,Ygrid]=meshgrid(linspace(L0,L0+L,NeleL+1),linspace(H0,H0+H,NeleH+1));
%Ygrid=Ygrid.*linspace(1,0.5,size(Ygrid,2));
nodeSeq=reshape(linspace(1,(NeleL+1)*(NeleH+1),(NeleL+1)*(NeleH+1)),NeleL+1,NeleH+1)';
eleSeq=reshape(linspace(1,(NeleL)*(NeleH),(NeleL)*(NeleH)),NeleL,NeleH)';
connect=zeros(NeleH*NeleL,4);
for eleH=1:NeleH
    for eleL=1:NeleL
    connect(eleSeq(eleH,eleL),:)=[nodeSeq(eleH,eleL),nodeSeq(eleH,eleL+1),nodeSeq(eleH+1,eleL+1),nodeSeq(eleH+1,eleL)];
    end
end
geometry.Econn=connect;
geometry.Coord=[reshape(Xgrid',(NeleL+1)*(NeleH+1),1),reshape(Ygrid',(NeleL+1)*(NeleH+1),1)];
geometry.Dof=transpose(reshape(1:(NeleL+1)*(NeleH+1)*2,2,(NeleL+1)*(NeleH+1)));
geometry.nodeSeq=nodeSeq;
geometry.eleSeq=eleSeq;
geometry.Xgrid=Xgrid;
geometry.Ygrid=Ygrid;
end

