function []=plotObj(quantity,geometry,varargin)
figure;
switch length(quantity)
    case geometry.NGPtotal
%MeshRefineFactor=varargin{1};
QuantityName=varargin{1};
%finemesh=refineGrid(geometry.Xgrid,geometry.Ygrid,MeshRefineFactor);
InterpolatQuantity=scatteredInterpolant(reshape(geometry.GPPhysiCoordx,geometry.NGPtotal,1),...
                        reshape(geometry.GPPhysiCoordy,geometry.NGPtotal,1),...
                        quantity(1:numel(geometry.Econn)));
% finemesh.Xgrid(102:end,[1:100,202:end])=nan;
% finemesh.Ygrid(102:end,[1:100,202:end])=nan;                    
surf(geometry.Xgrid,geometry.Ygrid,InterpolatQuantity(geometry.Xgrid,geometry.Ygrid))
    case size(geometry.Dof,1)
        InterpolatQuantity=scatteredInterpolant(geometry.Coord(:,1),...
                        geometry.Coord(:,2),...
                        quantity);
surf(geometry.Xgrid,geometry.Ygrid,InterpolatQuantity(geometry.Xgrid,geometry.Ygrid))
QuantityName=varargin{1};
end
   view(2)
      axis tight equal
      shading interp
      colormap(flipud(jet))
      colorbar
      title(QuantityName,'interpreter','latex')   
      xlabel('$x$[m]','interpreter','latex')
      ylabel('$y$[m]','interpreter','latex')
end