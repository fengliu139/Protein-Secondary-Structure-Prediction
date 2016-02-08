%convert protein sequences to numbers
for i=1:110
   seq(i).O=aa2int(seq(i).O);
end;
b=false;
%just to validate there are no special symbols in the sequences
for i=1:100
    if length(find(seq(i).O==25|seq(i).O==24|seq(i).O==21|seq(i).O==22|seq(i).O==23|seq(i).O==0))~=0
        b=true;
    end;
end;

for i=1:110
    seq(i).num_tran(1:3,1:3)=1;
    seq(i).num_emit(1:3,1:20)=1;
    last=state2num(seq(i).S(1));
    if seq(i).S(1)=='h'
            seq(i).num_emit(1,seq(i).O(1))=seq(i).num_emit(1,seq(i).O(1))+1;
    else if seq(i).S(1)=='e'
                seq(i).num_emit(2,seq(i).O(1))=seq(i).num_emit(2,seq(i).O(1))+1;
         else
                seq(i).num_emit(3,seq(i).O(1))=seq(i).num_emit(3,seq(i).O(1))+1;
         end;
    end;
    for j=2:length(seq(i).O)
        if seq(i).S(j)=='h'
            seq(i).num_emit(1,seq(i).O(j))=seq(i).num_emit(1,seq(i).O(j))+1;
        else if seq(i).S(j)=='e'
                seq(i).num_emit(2,seq(i).O(j))=seq(i).num_emit(2,seq(i).O(j))+1;
            else
                seq(i).num_emit(3,seq(i).O(j))=seq(i).num_emit(3,seq(i).O(j))+1;
            end;
        end;
        current=state2num(seq(i).S(j));
        seq(i).num_tran(last,current)=seq(i).num_tran(last,current)+1;
        last=current;
    end;
end;
 
%calculating the state transition probability matrix A and emission probability matrix B for every sequence
for i=1:110
    for k=1:3
        for I=1:3
                seq(i).A(k,I)=seq(i).num_tran(k,I)/sum(seq(i).num_tran(k,:));
        end;
        for b=1:20
                seq(i).B(k,b)=seq(i).num_emit(k,b)/sum(seq(i).num_emit(k,:));
         end;
    end;
end;

%calculating the combinational A and B
A(1:3,1:3)=0;
B(1:3,1:20)=0;
for i=1:110
    A=A+seq(i).A;
    B=B+seq(i).B;
end;
A=A/110;
B=B/110;


