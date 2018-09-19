function varargout=sortdata(order,varargin)
noi=length(varargin);
nod=length(order);
data=zeros(nod,noi+1);
orginialdirection=zeros(noi+1,1);
if size(order,2)==1;
    data(:,1)=order;
    orginialdirection(1)=1;
elseif size(order,1)==1;
    data(:,1)=order';
end

for i=1:noi
    if size(varargin{i},2)==1;
        data(:,i+1)=varargin{i};
        orginialdirection(i+1)=1;
    elseif size(varargin{i},1)==1;
        data(:,i+1)=varargin{i}';
    end
end
data=sortrows(data);
varargout=cell(noi+1,1);
for i=1:noi+1
    if orginialdirection(i)
        varargout{i}=data(:,i);
    else
        varargout{i}=data(:,i)';
    end
end
end