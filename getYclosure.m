%get Y with closure

function [Y]=getYclosure(strain,material,h)
array2tensor=@ (strain) [strain(1) strain(3); strain(3) strain(2)];
tensor2array=@ (tensor) [tensor(1), tensor(4), sqrt(2)*tensor(2)];
%tensor2voigt=@ (tensor) [tensor(1), tensor(4), 2*tensor(2)];
mu=material.mu;
lambda=material.lambda;
nGp=size(strain,1);
%strain=gpuArray(strain);
%strainTensor=zeros(2,2,size(strain,1));
Y=zeros(nGp,1);
for i=1:nGp
    strainTensor=array2tensor(strain(i,:));
    strainTr=strain(i,1)+strain(i,2);
    [Vec,Val]=eig(strainTensor);
    strainP=(Vec*max(Val,0)*transpose(Vec));
    strainN=strainTensor-strainP;
    
    
    strainP=tensor2array(strainP);
    strainN=tensor2array(strainN);
    
    YP=heaviside(strainTr)*(0.5*lambda*strainTr^2+mu*strainP*transpose(strainP));
    YN=h*(1-heaviside(strainTr))*(0.5*lambda*strainTr^2+mu*strainN*transpose(strainN));
    
    Y(i)=YP+YN;
    
end

% 
% newstrain=(max(strain,0)+min(strain,0).*h);
% 
% Y=0.5*sum(newstrain*C.*newstrain,2);
end


