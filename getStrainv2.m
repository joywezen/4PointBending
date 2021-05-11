function [straintransp]=getStrainv2(geometry,U)
straintransp=transpose(reshape(transpose(reshape(permute(geometry.B(:,:,:,1),[2,1,3]),8,12))*transpose(U(geometry.Edof(:,2:end))),3,geometry.NGPtotal));
end