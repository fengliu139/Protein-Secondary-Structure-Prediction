function seq=Parse(fn)
fd=fopen(fn);
i=0;
j=1;
line=fgetl(fd);
while ~feof(fd)
    if strcmp(line,'<>')
        i=i+1;
        line=fgetl(fd);
        while strcmp(line,'end')==0&strcmp(line,'<end>')==0&strcmp(line,'<>')==0
            
            seq(i).O(j)=line(1);
            seq(i).S(j)=line(3);
            j=j+1;
            line=fgetl(fd);
        end;
        if ~strcmp(line,'<>')
            line=fgetl(fd);
        end;
        j=1;
    else 
        line=fgetl(fd);
    end;
end;
fclose(fd);
for n=1:length(seq)
 seq(n).O=aa2int(seq(n).O);
end;