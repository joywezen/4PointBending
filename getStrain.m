function [strainbyEle,straintransp]=getStrain(geometry,U)
strainbyEle=zeros(3,4,size(geometry.Edof,1));
for i=1:size(geometry.Edof,1)
    strainbyEle(:,:,i)=reshape(transpose(reshape(permute(geometry.B(:,:,:,i),[2,1,3]),8,12))*U(geometry.Edof(i,2:end)),3,4);
end
%epsx=transpose(strain(1,:));
%epsy=transpose(strain(2,:));
%epsxy=transpose(strain(3,:));
straintransp=transpose(reshape((strainbyEle),3,i*4));
end