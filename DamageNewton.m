function [U_t1,D_t1,Y_t1,K_t1,Dflag,loop]=DamageNewton(material,geometry,F_ext,D_t0,U_t0,Y_t0,K_t0,tol,dt,randomness,sampleNr)
maxLoop=2e2;
rvsize=size(D_t0);
%rng('shuffle','twister');
switch randomness.NoiseType
    case 'prop-GS'
        rv1=randn(randomness.seeds{sampleNr},rvsize).*sqrt(dt)./3;
        rv2=0;
        rho=0;
    case 'prop-GL'
        rv1=randn(randomness.seeds{sampleNr},rvsize).*sqrt(dt)./3*10;
        rv2=0;
        rho=0;
    case 'prop-2GS'
        rv1=randn(randomness.seeds{sampleNr},rvsize).*sqrt(dt)./3*2;
        rv2=0;
        rho=0;
        
    case 'prop-0.5GS'
        rv1=randn(randomness.seeds{sampleNr},rvsize).*sqrt(dt)./3*0.5;
        rv2=0;
        rho=0;
    case 'non-prop-GS'
        
    case 'non-prop-GL'
    case 'prop-WBS'
    case 'prop-WBL'
    case 'noNoise'
        rv1=0;
        rv2=0;
        rho=0;
    otherwise
        fprintf('not implemented yet')
        return;
end
%if strcmp(NoiseType,'prop-GS')

damafun=@ (Y,DYDt,S,s1,s2) (Y./S).^s1.*(DYDt).^s2.*(1+(1-rho)*rv1./sqrt(dt))+(rho)*rv2./sqrt(dt);
%elseif strcmp(NoiseType,'prop-GL')
% rv1=randn(rvsize).*sqrt(dt)./3*10; %.*sqrt(dt)./3*10 is the maximum std
% rv2=0;%randn(rvsize).*sqrt(dt)./3;
% rho=0; %prop
%damafun=@ (DYDt,S,s) (DYDt./S).^s+(DYDt./S).^s*rv/sqrt(dt);
%damafun=@ (DYDt,S,s) (DYDt./S).^s;
%damafun=@ (DYDk,S,s) (DYD0k./S).^s;
% else
% damafun=@(Y,DYDt,S,s1,s2)(Y./S).^s1.*(DYDt).^s2;
% end
%damafun=@ (Y,DYDt,S,s1,s2) (Y./S).^s1.*(DYDt).^s2;
%% start Newton Method
U_k1=zeros(geometry.maxEdof,1);
r_k1=zeros(geometry.maxEdof,1);
r_t0(geometry.freedof,1)=K_t0(geometry.freedof,geometry.freedof)*U_t0(geometry.freedof)-F_ext(geometry.freedof);

%judge elastisity
delU_k1=K_t0(geometry.freedof,geometry.freedof)\(-r_t0(geometry.freedof,1));
U_k1(geometry.freedof)=U_t0(geometry.freedof)+delU_k1;

strain_k1=getStrainv2(geometry,U_k1);
% switch law
%     case 1
Y_k1=getYclosureGK(strain_k1,material); %#########has to be changed with getYclosure
%    case 2
%        Y_k1=getYclosure(strain_k1,material,h);
%end
if all(Y_k1(geometry.DamageGp)<=material.YD) && all(D_t0<=material.Dmax)  % all gauss points are inside ELASTIC region
    U_t1=U_k1;
    D_t1=D_t0;
    Dflag=0;
    loop=0;
    Y_t1=Y_k1;
    K_t1=K_t0;
elseif all(D_t0<=material.Dmax)                         % has somewhere DAMAGE evolution, determine bad locations
    badLocation=find(Y_k1>material.YD);
    %badLocation=intersect(badLocation,geometry.DamageGp);
    %start newton iteration
    loop=1;
    
    
    while 1
        DYDt=zeros(numel(geometry.Econn),1);
        DYDt(badLocation)=(Y_k1(badLocation)-max(repmat(material.YD,size(badLocation)),Y_t0(badLocation)))./dt;%0.5*sum(strain0*material.C.*(strain0),2);
        dD=damafun(max(Y_k1,0),max(DYDt,0),material.S,material.s1,material.s2).*dt;
        dD=max(dD,0);
        D_k1=D_t0+dD;
        %         switch law
        %             case 1
        %                 %tic
        K_k1=getTangentOperator(geometry,material,D_k1,strain_k1,'H');
        %toc
        %             case 2
        %
        %                 K_k1=getK(geometry,material,D_k1);
        %
        %         end
        r_k1(geometry.freedof,1)=K_k1(geometry.freedof,geometry.freedof)*U_k1(geometry.freedof)-F_ext(geometry.freedof);
        
        if all(D_k1<material.Dmax) && loop<=maxLoop
            
            if norm(r_k1(geometry.freedof))/norm(F_ext(geometry.freedof))<=tol && norm(delU_k1)/norm(U_k1(geometry.freedof))<=tol && all(D_k1>=D_t0)%any(abs(r(:)./F_ext(:))<tol)%any(D1-D0<tol) || any(abs(r)<tol)
                U_t1=U_k1;
                D_t1=D_k1;
                Dflag=1;
                Y_t1=Y_k1;
                K_t1=K_k1;
                return;
                %             elseif any(D_k1<D_t0)
                %                 U_t1=U_t0; %if accept reduced damage
                %                 D_t1=D_t0;
                %                 Dflag=2;
                %                 fprintf('error, damage reduced! BUT accepted \n')
                %                 end
                %                 report realistic damage
                %                 Y_t1=Y_t0;
                %                 K_t1=K_t0;
                %                 return;  % break the while loop and accept D_t1
            else  %continue next while loop
                delU_k1=K_k1(geometry.freedof,geometry.freedof)\(-r_k1(geometry.freedof,1));
                U_k1(geometry.freedof)=U_k1(geometry.freedof)+delU_k1;
                strain_k1=getStrainv2(geometry,U_k1);
                %                 switch law
                %                     case 1
                Y_k1=getYclosureGK(strain_k1,material);%0.5*sum(strain_k1*material.C.*(strain_k1),2);
                %                     case 2
                %                         Y_k1=getYclosure(strain_k1,material,h);
                %                 end
                
                %Y_k1=getYclosureGK(strain_k1,material);%0.5*sum(strain_k1*material.C.*(strain_k1),2);
                loop=loop+1;
            end
            
            
            
        elseif loop>maxLoop && all(D_k1<material.Dmax)
            fprintf('the number of newton iteration is greater than material.Dmax %3.0f! accept value from last iteration \n',maxLoop)
            U_t1=U_k1;
            D_t1=D_k1;
            Dflag=3;
            Y_t1=Y_k1;
            K_t1=K_k1;
            return;
            
        else % any(D_k1>=1)
            fprintf('damage is greater than material.Dmax %3.2f ! keep previous value \n', material.Dmax)
            D_k1(D_k1>material.Dmax)=material.Dmax;
            U_t1=U_k1;
            D_t1=D_k1;
            Y_t1=Y_k1;
            K_t1=K_k1;
            Dflag=4;
            return;
            
        end
        %
    end
else
    fprintf('previous value is greater than %3.2f, newton loop will not contitinue, keep previous value \n',material.Dmax)
    U_t1=U_t0;
    D_t1=D_t0;
    Dflag=1;
    loop=0;
end

end
