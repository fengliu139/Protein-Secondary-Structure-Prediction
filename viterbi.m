function [test,v]=viterbi(test,a,e,pie)

for n=1:17
%set the first column of v and ptr matrix
x=test(n).O(1);
for i=1:3
    v(i,1)=log2(e(i,x))+1+log2(pie(i));
    ptr(i,1)=0;
end;

%iteration of viterbi
l=length(test(n).O);
for i=2:l
x=test(n).O(i);
    for j=1:3
        for k=1:3
            temp(k)=v(k,i-1)+log2(a(k,j));
        end;
        [m,kk]=max(temp);
        v(j,i)=m+log2(e(j,x));
        ptr(j,i)=kk;
    end;
end;

%trace back
[m,p] = max(v(:,l));
test(n).s(1:l)=0;
test(n).s(l)=p;
for i=(l-1):-1:1
    test(n).s(i)=ptr(test(n).s(i+1),i+1);
end;
for i=1:l
    switch test(n).s(i)
        case 1
            test(n).ss(i)='h';
        case 2
            test(n).ss(i)='e';
        otherwise
            test(n).ss(i)='_';
    end;
end;

end;