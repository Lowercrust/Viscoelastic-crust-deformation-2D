% calculate 

function F=FPT(p,T,mesh,type)
% input mesh with 2 face geometry
% input argument p(pressure[Pa]) T(temperature[K])
% output argument F [MPa^{n}*s]
% global rhe rhe1 R L 
global modc fi
crhe=modc.rhe{1};
mrhe=modc.rhe{2};
Adis_ano=10^crhe.logA; % anorthite dislocation
Adif_ano=10^crhe.logAd; % anorthite diffusion
Adis_oli=10^mrhe.logA;% olivine dislocation
Adif_oli=10^mrhe.logAd;% olivine diffusion
% interp P and T to the center of triangle
% [p1,~,t1]=meshToPet(mesh);
p1=mesh.Nodes;
% c = Composite(2);
% c{1}=p';
% c{2}=T';
% spmd(2)
%     c3=pdeintrpgpu(p1,t1,c);
% end
% p=c3{1};
% T=c3{2};

% p=pdeintrpgpu(p1,t1,p');
% T=pdeintrpgpu(p1,t1,T');
p1i=fi{1};
p2i=fi{2};
% center coordinate of triangles
%[x,y]=tricenter(p1);
% face1tri=find(t1(4,:)==1); % crust
% face2tri=find(t1(4,:)==2); % mantle
F=zeros(1,length(p1));
% dislocation creep
if strcmp(type,'disl') && crhe.m==0 && mrhe.m==0
    F(1,p1i)=(1./(Adis_ano.*(fugacity(p(p1i),T(p1i))).^crhe.r)).*exp((crhe.Q+p(p1i).*crhe.V)./(modc.R.*T(p1i)));
    F(1,p2i)=(1./(Adis_oli.*(fugacity(p(p2i),T(p2i))).^mrhe.r)).*exp((mrhe.Q+p(p2i).*mrhe.V)./(modc.R.*T(p2i)));
elseif strcmp(type,'diff') && crhe.md~=0 && mrhe.md~=0
    %only for the situation when grain size is known(e.g. assuming
    %constant grain size)
    L=modc.L;
    F=(L.^rhe.md/(Adif_ano.*(fugacity(p,T)).^rhe.rd)).*exp((rhe.Qd+p.*rhe.Vd)./(R.*T));
else
    error('input disl or diff for type (FPT(p,T,type))')
end
% F=pdeprtnigpu(p1,t1,F)';
end