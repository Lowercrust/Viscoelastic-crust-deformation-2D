function uinit = initaltemperature(locations)
global inittempf
M = length(locations.x);
%uinit = zeros(1,M);
uinit=inittempf(locations.x,locations.y);