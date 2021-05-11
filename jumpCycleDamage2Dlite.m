function [Dout,N_f,DCycleBeginOut,DCycleEndOut,jumpCycle,jumpAcc]=jumpCycleDamage2Dlite(material,geometry,maxCycle,Kappa,Ndt,dt,F_ext,randomness,sampleNr)
%initialize
cycle=1;
tol=1e-6;
DCycleBegin=zeros(geometry.NGPtotal,1);
DCycleEndOut=zeros(geometry.NGPtotal,Kappa*2);
DCycleBeginOut=zeros(geometry.NGPtotal,Kappa*2);
D_t0=DCycleBegin;
jumpCycle=zeros(1,Kappa*2);
%Dout=zeros(geometry.NGPtotal,Kappa*2);
U_t0=zeros(geometry.maxEdof,1);
Y_t0=zeros(geometry.NGPtotal,1);
K_t0=getK(geometry,material,D_t0);
jumpAcc=1;
startT=tic;
while 1
    %% loading phase
    for i=(2:Ndt/2+1)+(cycle-1)*Ndt
        %fprintf('step %d\n',i-1)
        [U,D,Y,K,Dflag,~]=DamageNewton(material,geometry,F_ext(:,i-(cycle-1)*Ndt),D_t0,U_t0,Y_t0,K_t0,tol,dt,randomness,sampleNr);
        if Dflag==4 % D already over Dmax
         %   fprintf('Warning! damage evaluation is terminated, keep the last effective value \n')
            Dout=D;
            N_f=cycle;
            return;
        end
        %next time step
        D_t0=D;
        U_t0=U;
        K_t0=K;
        Y_t0=Y;
        %Dmax(i)=max(D);
    end
    
    %% unloading phase
    U_t0=zeros(geometry.maxEdof,1);
    Y_t0=zeros(geometry.NGPtotal,1);
    DCycleEnd=D;
    %% output CycleEnd and CycleBegin
    DCycleEndOut(:,jumpAcc)=DCycleEnd;
    DCycleBeginOut(:,jumpAcc)=DCycleBegin;
    %% jump phase
    dDdN=DCycleEnd-DCycleBegin;
    jumpCycle(jumpAcc)=floor(min((material.Dmax/Kappa)./dDdN));
    jumpEndTime=toc(startT);
    if jumpCycle(jumpAcc)==inf 
        %fprintf('Warning! No damage, increase the load \n')
        Dout=D;
        N_f=cycle;
        return;
    elseif jumpCycle(jumpAcc)>=maxCycle
        Dout=D;
        N_f=jumpCycle(jumpAcc)+cycle;
        return;
    elseif jumpCycle(jumpAcc)+cycle>=maxCycle
        jumpCycle(jumpAcc)=maxCycle-cycle;
        DCycleEnd=DCycleEnd+dDdN*jumpCycle(jumpAcc);
        %judge with Dmax
        if max(DCycleEnd)<material.Dmax
            Dout=DcycleEnd;
            N_f=maxCycle; %finish maxCycle but model still not fail
        else
            JumpCycleCorrection=-ceil(min((max(DCycleEnd)-material.Dmax)./dDdN));
            Dout=DCycleEnd+dDdN*JumpCycleCorrection;
            N_f=cycle+JumpCycleCorrection; %within maxCycle, model fail
            
        end
        return;
    else
        DCycleEnd=DCycleEnd+dDdN*jumpCycle(jumpAcc);
        %% attention debug on 17.April %%
         cycle=cycle+jumpCycle(jumpAcc)+1; % should be account firt before judge
         jumpAcc=jumpAcc+1;%% should be account firt before judge
        %%
        
        if max(DCycleEnd)<material.Dmax
            N_f=0; %within maxCycle, model not fail, update and continue full computation
            DCycleBegin=DCycleEnd;
            D_t0=DCycleBegin;
            Dout=DCycleEnd;
             fprintf('%0.1f[s] %0.1f%% elapsed, %0.1f[s] till finish \n',jumpEndTime,jumpAcc/Kappa*100,jumpEndTime/(jumpAcc/Kappa)-jumpEndTime)
             %fprintf('%d of %d cycles elapsed, the following %d cycles will be avoided \n',cycle,maxCycle,jumpCycle(jumpAcc));   
            
        else
            JumpCycleCorrection=-ceil(min((max(DCycleEnd)-material.Dmax)./dDdN));
            Dout=DCycleEnd+dDdN*JumpCycleCorrection;
            N_f=cycle+JumpCycleCorrection;
            jumpCycle(jumpAcc-1)=jumpCycle(jumpAcc-1)+JumpCycleCorrection;
            return;
        end
    end
end
end
