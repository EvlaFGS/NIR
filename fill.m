1;

function t = account_for_periodic (x,ed)
t=x;
if (x<-ed)
t=x+2*ed;
elseif (x>ed)
t=x-2*ed;
endif
end

xx=[];
yy=[];
zz=[];
ii=[];

 A = dlmread('Coord.txt','\t',1,0)

prompt = "Enter A_C in mkm: ";
A_C = input(prompt)
prompt = "Enter B in mkm: ";
B = input(prompt)
prompt = "Enter edge of a cube: ";
EDGE=input(prompt)

ed=EDGE/2.0;
for i=1:length(A(:,1))
for x=-ed-B:0.5:ed+B
for y=-ed-B:0.5:ed+B
for z=-ed-B:0.5:ed+B
if ((x-A(i,1)).^2/(A_C.^2)+(y-A(i,2)).^2/(B.^2)+(z-A(i,3)).^2/(A_C.^2)<=1.0)
x=account_for_periodic(x,ed);
y=account_for_periodic(y,ed);
z=account_for_periodic(z,ed);
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
