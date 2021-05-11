function [Y] = getYclosureGK(strain,material)
% switch type
%     case 'fullY'
% %array2tensor=@ (strain) [strain(1) strain(3); strain(3) strain(2)];
% %tensor2array=@ (tensor) [tensor(1), tensor(4), tensor(2)];
% getStrainDeviator = @ (Strain,StrainTr) [[Strain(:,1) Strain(:,2)]-1/3*StrainTr, Strain(:,3)];
% StrainDeviator2Squar=@ (StrainDeviatorArray) StrainDeviatorArray(:,1).^2+StrainDeviatorArray(:,2).^2+2*StrainDeviatorArray(:,3).^2;
% G=material.G;
% K=material.K;
% 
% 
%     StrainTr=strain(:,1)+strain(:,2);
%     StrainDeviatorArray=getStrainDeviator(strain,StrainTr);
%     Y=G*StrainDeviator2Squar(StrainDeviatorArray)+K/2.*(max(StrainTr,0)).^2;
%     %Y=K/2.*(max(StrainTr,0)).^2;
%     case 'H'
    K=material.K;
    StrainTr=strain(:,1)+strain(:,2);
    Y=K/2.*(max(StrainTr,0)).^2;
% end
end