function imwritematrix(A,colormap,datalimit)
    maxA=max(max(A));
    minA=min(min(A));
    colorlength=length(colormap);
    % diverging colormap
    A(A<datalimit(1))=datalimit(1);
    A(A>datalimit(2))=datalimit(2);
    if datalimit(1)<0 && datalimit(2)>0
        rangeplus=datalimit(2);
        rangeminus=abs(datalimit(1));
        A(A>0)=A(A>0)*colorlength/(2*rangeplus);
        A(A<0)=A(A<0)*colorlength/(2*rangeminus);
        A=A+colorlength/2;
    else % sequential colormap
        range=abs(datalimit(2)-datalimit(1));
        A=(A-datalimit(1))*colorlength/range;
    end
    imshow(A,colormap);
end