function [NoDamageGp] = findeNearGp(geometry,NoDamageNode)
iGp=0;
for i=1:size(geometry.Econn,1)
if any(ismember(NoDamageNode,geometry.Econn(i,:)))
    iGp=iGp+1;
    NoDamageGp(iGp,:)=[(i-1)*4+1:(i-1)*4+4];
end
end
NoDamageGp=reshape(NoDamageGp',numel(NoDamageGp),1);
end