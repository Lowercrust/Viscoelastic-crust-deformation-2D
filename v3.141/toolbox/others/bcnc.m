function bcMatrix = bcnc(problem,region,state)
global vbc mod
if vbc(region.y)<0
    bcMatrix=0;
elseif vbc(region.y)*(365*24*3600*1000)>mod.v0
    bcMatrix=mod.v0/(365*24*3600*1000);
else
    bcMatrix=vbc(region.y);
end
end