xx=[];
yy=[];
zz=[];
ii=[];

 A = dlmread('Coord.txt','\t',1,0)


for i=1:1:30
for x=-10:0.6:10
for y=-10:0.6:10
for z=-10:0.6:10
if ((x-A(i,1)).^2/4.0+(y-A(i,2)).^2/4.0+(z-A(i,3)).^2/4.0<1.0)
xx=[xx x];
yy=[yy y];
zz=[zz z];
ii=[ii i];
end

end
end
end
end

scatter3(xx,yy,zz,[],ii)
pause
