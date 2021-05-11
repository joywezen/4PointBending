function [fine] = refineGrid(Xgrid, Ygrid, varargin)
if ~isempty(varargin)
[fine.Xgrid,fine.Ygrid]=meshgrid(linspace(Xgrid(1),Xgrid(end),(size(Xgrid,2)-1)*varargin{:}+1),...
    linspace(Ygrid(1),Ygrid(end),(size(Ygrid,1)-1)*varargin{:}+1));
else
    fine.Xgrid=Xgrid;
    fine.Ygrid=Ygrid;
end
end