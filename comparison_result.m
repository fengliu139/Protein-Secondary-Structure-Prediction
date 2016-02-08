for n=1:17
k=(test(n).ss==test(n).S);
result(n)=length(find(k==1))/length(test(n).S);
end;
sum(result)/17