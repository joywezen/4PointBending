function [vonMises,sigx,sigy,tauxy]=getvonMises(geometry,material,D,U)
stress=zeros(3,4,size(geometry.Edof,1));
for j=1:size(geometry.Edof,1)
stress(:,:,j)=reshape(blkdiag(hooke(1,material.E,material.nv)*(1-D((j-1)*4+1)),hooke(1,material.E,material.nv)*(1-D((j-1)*4+2)),...
    hooke(1,material.E,material.nv)*(1-D((j-1)*4+3)),hooke(1,material.E,material.nv)*(1-D((j-1)*4+4)))*...
    transpose(reshape(permute(geometry.B(:,:,:,j),[2,1,3]),8,12))*U(geometry.Edof(j,2:end)),3,4);
end
sigx=transpose(stress(1,:));
sigy=transpose(stress(2,:));
tauxy=transpose(stress(3,:));
%vonMises=reshape((sigx.^2+sigy.^2-sigx.*sigy+3*tauxy.^2).^(0.5), 4,size(geometry.Edof,1));
vonMises=((sigx.^2+sigy.^2-sigx.*sigy+3*tauxy.^2).^(0.5));
end