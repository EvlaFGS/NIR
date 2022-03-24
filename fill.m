xx=[];
yy=[];
zz=[];
ii=[];

 A = dlmread('Coord.txt','\t',1,0)

prompt = "Enter A_C in mkm: ";
A_C = input(prompt)
prompt = "Enter B in mkm: ";
B = input(prompt)

ed=18.155/2.0;
for i=1:1:30
for x=-ed:0.5:ed
for y=-ed:0.5:ed
for z=-ed:0.5:ed
if ((x-A(i,1)).^2/(A_C.^2)+(y-A(i,2)).^2/(B.^2)+(z-A(i,3)).^2/(A_C.^2)<1.0)
xx=[xx x];
yy=[yy y];
zz=[zz z];
ii=[ii i];
end

end
end
end
end

pp=[xx; yy; zz]';
dlmwrite("voxels_octave.txt",pp,' ');

scatter3(xx,yy,zz,[],ii)
pause
