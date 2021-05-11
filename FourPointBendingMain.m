%----- System matrices -----
clear all
close all
clc
%s = RandStream('mt19937ar','Seed','shuffle');
%RandStream.setGlobalStream(s);
%addpath(genpath('../../external_functions'));
%addpath('../../OneDBarFEM');
%addpath(genpath('../'));
%environment_setting
material.S=1.510e3;%1868;%2.245e3%1.731e3%7.438e2 %between 1.868e3 and 1.867e3
material.s1=11.05;%12.05
material.s2=1;
material.fr=6.05e6;
material.sigmaD=material.fr*0.1;
material.E=42e9;
material.nv=0.1;
material.G=material.E/2/(1+material.nv);
material.K=material.E/3/(1-2*material.nv);
material.mu=material.G;
material.lambda=material.K-(2*material.G/3);
material.C=hooke(1,material.E,material.nv);
material.YD=material.sigmaD^2/2/material.E;



%Ep=[1, t]; %1 for plane stress, 2 for plane strain
%disp('PRESS ENTER TO CONTINUE'); pause; clf; 


%%3 for 2d problem,number of freedoms per element,number of gp per element,number of element;

% B=zeros(3,8,4,2);
% for i=1:size(geometry.Econn,1)
% [Ke,~,B(:,:,:,i)]=planre(geometry.Ex(i,:),geometry.Ey(i,:),Ep,material.C,[],D(:,i));
% K=assem(geometry.Edof(i,:),K,Ke);
% end
geometry.ProblemDimension=2;
geometry.H=1e-1;
geometry.L=5e-1;
geometry.L0=0;
geometry.H0=0;
geometry.NeleH=8;
geometry.NeleL=40;
geometry.Nele=geometry.NeleH*geometry.NeleL;
geometry.directions=[1,2];
geometry.fixCoord=[0.025,0;0.475,0];%[0.18,0;0.2,0;0.8,0;0.82,0];[zeros(11,1),linspace(0,0.1,11)'];%[
geometry.loadCoord=[0.2,0.1;0.3,0.1];%[0.125,0.1;0.375,0.1];%[0.4,0.1;0.42,0.1;0.58,0.1;0.6,0.1];
[geometry]=setObject(2,geometry);
%close all
[geometry]=getBmat(geometry);
[geometry.GPPhysiCoordx,geometry.GPPhysiCoordy]=getPhysiCoord(geometry);
F=@(sigma, b,d,L,Li) sigma*2*b*d^2/(L-Li)/3;

Fmax=F(material.fr,geometry.H,geometry.t,diff(geometry.fixCoord(:,1)),diff(geometry.loadCoord(:,1)));%1e-3*material.fr;
Fmin=0.1*Fmax;
f=10;
Ndt=400; % 100 time steps per cycle
maxCycle=2e8;
t_vec=linspace(0,1/f,Ndt+1);
dt=t_vec(2)-t_vec(1);
%th_A=Fmin;
%loadtype='saw';

force_genSin=@(Fmax,Fmin,f,t) -0.5*(Fmax-Fmin)*cos(2*pi*f*t)+(0.5*(Fmax-Fmin)+Fmin);
%force_genTri=@(Fmax,Fmin,f,t) sawtooth(t*f*2*pi,0.5)*(Fmax-Fmin)/2+(Fmin+(Fmax-Fmin)/2); %%triangle wave


F_ext_Factor=force_genSin(Fmax,Fmin,f,t_vec);%force_generator(Fmax,Fmin,f,th_A,loadtype,t_vec);

F_ext=zeros(geometry.maxEdof,length(t_vec));
F_ext(geometry.Dof(geometry.loadnodes,geometry.directions(1)),1:length(t_vec))=0;
F_ext(geometry.Dof(geometry.loadnodes,geometry.directions(2)),1:length(t_vec))=-repmat(F_ext_Factor,length(geometry.loadnodes),1);%-0.1*material.fr*0.005;


loadLevelFactor=[0.7];
loadLevels=length(loadLevelFactor);
material.Dmax=0.3;
Kappa=100;
%jumpCycle=0;
%cycle=1;

%random noises parameters
mu = 1;
sig2 =1;
fun = @(k)sig2/mu^2-gamma(1+2./k)./gamma(1+1./k).^2+1;
k0 = 1;               % Initial guess
material.k = fzero(fun,k0);       % Solve for k
material.lam = mu/gamma(1+1/material.k); % Substitue to find lambda
numcore=feature('numcores');
switch numcore
    case 6
        Nsamples=150;
    case 8
Nsamples=200;
    case 16
        Nsamples=100;
    case 4
        Nsamples=100;
end

DEnd=zeros(geometry.NGPtotal,Nsamples,loadLevels);
DCycleBegin=zeros(geometry.NGPtotal,Kappa*2,Nsamples,loadLevels);
DCycleEnd=zeros(geometry.NGPtotal,Kappa*2,Nsamples,loadLevels);
N_f=zeros(Nsamples,loadLevels);
jumpAcc=zeros(Nsamples,loadLevels);
jumpCycle=zeros(Kappa*2,Nsamples,loadLevels);
randomness.NoiseType='prop-2GS';
randomness.seeds=RandStream.create('mrg32k3a','NumStreams',Nsamples,'Seed','shuffle','CellOutput',true);
tic
%
formatOut = 'dd.mm.yyyy_HH:MM';
[~,hostname]=unix('hostname');
tempDateTime=datestr(datetime,formatOut);
for Ls=1:loadLevels
    F_extLevelS=F_ext.*loadLevelFactor(Ls);
    fprintf('**load level= %0.2f** \n',loadLevelFactor(Ls))
    
    %parfor_progress(Nsamples);
    parfor sampleNr=1:Nsamples    
       
        [DEnd(:,sampleNr,Ls),N_f(sampleNr,Ls),DCycleBegin(:,:,sampleNr,Ls),...
            DCycleEnd(:,:,sampleNr,Ls),jumpCycle(:,sampleNr,Ls),jumpAcc(sampleNr,Ls)]...
            = jumpCycleDamage2Dlite(material,geometry,maxCycle,Kappa,Ndt,dt,F_extLevelS,randomness,sampleNr);
        
    %    parfor_progress;
    end
   % parfor_progress(0);

    

%save (strcat(pwd,'/tempResultsMat/',hostname(1:6),'_',num2str(Nsamples),'Samples_','LoadLevel',num2str(loadLevelFactor(Ls)),NoiseType,tempDateTime,'.mat'),'-v7.3')

end


 ElapsT=toc;
% 
% 
% 
 beep
 save (strcat(pwd,'/resultsMat/',hostname(1:6),'_',num2str(Nsamples),randomness.NoiseType,tempDateTime,'.mat'),'-v7.3')
 exit

