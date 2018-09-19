function out=removeinf(in)
in(in==inf)=0;
in(in==0)=max(in);
in(in==-inf)=0;
in(in==0)=min(in);
out=in;
end