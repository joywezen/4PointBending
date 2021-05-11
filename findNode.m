function [Node]= findNode(geometry,targetCoord)
NoN=size(targetCoord,1);
Node=zeros(NoN,1);
for i=1:NoN
Node(i)=intersect(find(geometry.Coord(:,1)==targetCoord(i,1)),find(geometry.Coord(:,2)==targetCoord(i,2)));
end
end