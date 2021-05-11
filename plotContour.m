function []=plotContour(geometry,material,D,U,varargin)
environment_setting
formatOut = 'mm-dd-yy-HH-MM-SS-AM';
if isempty(varargin)
    fprintf('error, define target to be plotted, damage, damageUncertainty, displacement, stress, etc.\n')
    return;
else
    [vonMises,sigx,sigy,tauxy]=(getvonMises(geometry,material,D,U));
    
    
        finemesh=refineGrid(geometry.Xgrid,geometry.Ygrid,varargin{2});
    
    if strcmp( varargin{1}, 'damage')
      figure;
      Fdama=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,numel(geometry.GPPhysiCoordx),1),...
                        reshape(geometry.GPPhysiCoordy,numel(geometry.GPPhysiCoordy),1),...
                        D(1:numel(geometry.Econn)),'natural');
                    
      surf(finemesh.Xgrid,finemesh.Ygrid,Fdama(finemesh.Xgrid,finemesh.Ygrid)) 
      hold on
      
      view(2)
      axis tight equal
      colormap jet
      shading interp
      colorbar
      %hbad=plot3(geometry.GPPhysiCoordx(badLocation),geometry.GPPhysiCoordy(badLocation),[1,1],'o','color',colorsetgreen(6,:),'MarkerSize',10);
      caxis([0 material.Dmax*0.6])
      %pause(0.01)
      title('damage')
      legend off
      if strcmp(varargin{3},'tikz')
          filename=datestr(datetime,formatOut);
          filename=filename(~isspace(filename));
      convertTikz(tikzpath,strcat('2DbendingFinalDamageRandom',filename,'.tikz'),100)
      end
    elseif strcmp(varargin{1},'damageMean')
        figure;
        Dmea=mean(D,2); % size(D)=number of gauss points * number of samples
        FDstd=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,numel(geometry.GPPhysiCoordx),1),...
                        reshape(geometry.GPPhysiCoordy,numel(geometry.GPPhysiCoordy),1),...
                        Dmea(1:numel(geometry.Econn)));
                    
        surf(finemesh.Xgrid,finemesh.Ygrid,FDstd(finemesh.Xgrid,finemesh.Ygrid))
         view(2)
      axis tight equal
      shading interp
      colormap jet
      colorbar
      title('mean of damage')
      %caxis([0 0.3])
      if strcmp(varargin{3},'tikz')
          filename=datestr(datetime,formatOut);
          filename=filename(~isspace(filename));
      convertTikz(tikzpath,strcat('2DbendingDamageMean',filename,'.tikz'),100)   
      end
    elseif strcmp(varargin{1}, 'damageUncertainty')
        figure;
        Dstd=std(D,0,2); % size(D)=number of gauss points * number of samples
        FDstd=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,numel(geometry.GPPhysiCoordx),1),...
                        reshape(geometry.GPPhysiCoordy,numel(geometry.GPPhysiCoordy),1),...
                        Dstd(1:numel(geometry.Econn)));
                    
        surf(finemesh.Xgrid,finemesh.Ygrid,FDstd(finemesh.Xgrid,finemesh.Ygrid))
         view(2)
      axis tight equal
      shading interp
      colormap summer
      colorbar
      title('damage uncertainty')
      caxis([0 0.0062])
      if strcmp(varargin{3},'tikz')
          filename=datestr(datetime,formatOut);
          filename=filename(~isspace(filename));
      convertTikz(tikzpath,strcat('2DbendingDamageUncertainty',filename,'.tikz'),100)   
      end
      
    elseif strcmp(varargin{1}, 'Mises')
        figure;
        FvonMises=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,numel(geometry.GPPhysiCoordx),1),...
                        reshape(geometry.GPPhysiCoordy,numel(geometry.GPPhysiCoordy),1),...
                        vonMises(1:numel(geometry.Econn))); 
        contourf(finemesh.Xgrid,finemesh.Ygrid,FvonMises(finemesh.Xgrid,finemesh.Ygrid))  
        view(2)
        axis equal
        shading interp
        colorbar
        title('von Mises')            
                    
    elseif strcmp(varargin{1},'stressX')
        figure;
        Fsigx=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,numel(geometry.GPPhysiCoordx),1),...
                        reshape(geometry.GPPhysiCoordy,numel(geometry.GPPhysiCoordy),1),...
                        sigx(1:numel(geometry.Econn)));
         contourf(finemesh.Xgrid,finemesh.Ygrid,Fsigx(finemesh.Xgrid,finemesh.Ygrid))  
         view(2)
         axis equal
         shading interp
         colorbar
         title('$\sigma_x$','interpreter','latex')
         
    elseif strcmp(varargin{1},'stressY')
        figure;
        Fsigy=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,numel(geometry.GPPhysiCoordx),1),...
                        reshape(geometry.GPPhysiCoordy,numel(geometry.GPPhysiCoordy),1),...
                        sigy(1:numel(geometry.Econn)));
                    
         contourf(finemesh.Xgrid,finemesh.Ygrid,Fsigy(finemesh.Xgrid,finemesh.Ygrid))  
         view(2)
         axis tight equal       
         shading interp
         colorbar
         title('$\sigma_y$','interpreter','latex')
         
    elseif strcmp(varargin{1},'stressXY') 
        figure;
        Ftauxy=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,numel(geometry.GPPhysiCoordx),1),...
                        reshape(geometry.GPPhysiCoordy,numel(geometry.GPPhysiCoordy),1),...
                        tauxy(1:numel(geometry.Econn)));   
         contourf(finemesh.Xgrid,finemesh.Ygrid,Ftauxy(finemesh.Xgrid,finemesh.Ygrid))  
         view(2)
         axis tight equal
         shading interp
         colorbar
         title('$\tau_{xy}$','interpreter','latex')
    else 
        fprintf('error, define target to be plotted, damage, damageUncertainty, displacement, stress, etc. \n ')
        return;
    end
                    
     
    
     
     
                   
 %subplot(2,2,1)              
 
  %subplot(2,2,2)              
 
 %subplot(2,2,3)              

  %subplot(2,2,4)              
 
% pause(0.01)
% set(gcf, 'Position', get(0, 'Screensize'));

end
