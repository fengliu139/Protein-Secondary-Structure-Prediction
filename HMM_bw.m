%set initial parameters for the HMM
 the following is just an example
for n=1:110
    seq(n).A=[0.95 0.05 0.05;0.05 0.95 0.05; 0.05 0.05 0.95];
    seq(n).B(1:3,1:20)=0.05;
    seq(n).pie=[0.01 0.03 0.96];
end;

%calculate the length of all the training sequences
for n=1:110
    l(n)=length(seq(n).O);
end;

%the first iteration of BW for every sequences. It actually can be merged to the following iterations
%it was used to validate the BW implementation was corrent
for n=1:110
    for i=1:3
        seq(n).fv(1,i)=seq(n).pie(i)*seq(n).B(i,seq(n).O(1));
        seq(n).bv(l(n),i)=1;
    end;
    seq(n).Z(1)=sum(seq(n).fv(1,:));
    seq(n).fv(1,:)=seq(n).fv(1,:)/seq(n).Z(1);
    
    for t=1:(l(n)-1)
        for j=1:3
            seq(n).fv(t+1,j)=(seq(n).fv(t,:)*seq(n).A(:,j))*seq(n).B(j,seq(n).O(t+1));
            
        end;
        seq(n).Z(t+1)=sum(seq(n).fv(t+1,:));
        seq(n).fv(t+1,:)=seq(n).fv(t+1,:)/seq(n).Z(t+1);      
        
    end;
    for t=(l(n)-1):-1:1
        for i=1:3
            seq(n).bv(t,i)=seq(n).A(i,1)*seq(n).B(1,seq(n).O(t+1))*seq(n).bv(t+1,1)+seq(n).A(i,2)*seq(n).B(2,seq(n).O(t+1))*seq(n).bv(t+1,2)+seq(n).A(i,3)*seq(n).B(3,seq(n).O(t+1))*seq(n).bv(t+1,3);
        end;
        seq(n).bv(t,:)=seq(n).bv(t,:)/seq(n).Z(t+1);
        
    end;
    
    for t=1:l(n)
        theta(n,t)=seq(n).fv(t,1)*seq(n).bv(t,1)+seq(n).fv(t,2)*seq(n).bv(t,2)+seq(n).fv(t,3)*seq(n).bv(t,3);
    end;
    
    
    for t=1:l(n)-1
        for i=1:3
            for j=1:3
                seq(n).E(t,i,j)=(seq(n).fv(t,i)*seq(n).A(i,j)*seq(n).B(j,seq(n).O(t+1))*seq(n).bv(t+1,j))/(seq(n).fv(t,1)*seq(n).bv(t,1)+seq(n).fv(t,2)*seq(n).bv(t,2)+seq(n).fv(t,3)*seq(n).bv(t,3));
                seq(n).E(t,i,j)=seq(n).E(t,i,j)/seq(n).Z(t+1);
                
            end;
            seq(n).R(t,i)=seq(n).E(t,i,1)+seq(n).E(t,i,2)+seq(n).E(t,i,3);
        end;
        
    end;
    for i=1:3
        seq(n).R(l(n),i)=seq(n).fv(l(n),i)/sum(seq(n).fv(l(n),:));
    end;
    
    
    for i=1:3
        sumR=sum(seq(n).R(1:l(n)-1,i));
        for j=1:3
            seq(n).value(i,j)=sum(seq(n).E(:,i,j))/sumR;
            
        end;
    end;
    
    
    
    for i=1:3
        sumR=sum(seq(n).R(:,i));
        
        for j=1:20
            
            seq(n).value2(i,j)=sum(seq(n).R(seq(n).O==j,i))/sumR;
        end;
        
    end;
    
    seq(n).pie2=seq(n).R(1,:);
    
    seq(n).current=0;
    for i=1:l(n)
        seq(n).current=seq(n).current+log10(seq(n).Z(i));
    end;
end;


%BW iterations. 
for n=1:110
    seq(n).last=10000;
    %For every sequence, iterating until last P(O|λ) is close enough to the current P(O|λ)
    while (abs(seq(n).last-seq(n).current)>0.0001)
        seq(n).A=seq(n).value;
        seq(n).B=seq(n).value2;
        for i=1:3
            seq(n).fv(1,i)=seq(n).pie(i)*seq(n).B(i,seq(n).O(1));
            seq(n).bv(l(n),i)=1;
        end;
	 
        seq(n).Z(1)=sum(seq(n).fv(1,:));  %seq(n).Z is used to normalize fv and bv so they won't underflow
        seq(n).fv(1,:)=seq(n).fv(1,:)/seq(n).Z(1);
        
        for t=1:(l(n)-1)
            for j=1:3
                seq(n).fv(t+1,j)=(seq(n).fv(t,:)*seq(n).A(:,j))*seq(n).B(j,seq(n).O(t+1));
                
            end;
            seq(n).Z(t+1)=sum(seq(n).fv(t+1,:));
            seq(n).fv(t+1,:)=seq(n).fv(t+1,:)/seq(n).Z(t+1);
            
            
        end;
        for t=(l(n)-1):-1:1
            for i=1:3
                seq(n).bv(t,i)=seq(n).A(i,1)*seq(n).B(1,seq(n).O(t+1))*seq(n).bv(t+1,1)+seq(n).A(i,2)*seq(n).B(2,seq(n).O(t+1))*seq(n).bv(t+1,2)+seq(n).A(i,3)*seq(n).B(3,seq(n).O(t+1))*seq(n).bv(t+1,3);
            end;
            seq(n).bv(t,:)=seq(n).bv(t,:)/seq(n).Z(t+1);
            
        end;
        
        for t=1:l(n)
            theta(n,t)=seq(n).fv(t,1)*seq(n).bv(t,1)+seq(n).fv(t,2)*seq(n).bv(t,2)+seq(n).fv(t,3)*seq(n).bv(t,3);
        end;
        
        
        for t=1:l(n)-1
            for i=1:3
                for j=1:3
                    seq(n).E(t,i,j)=(seq(n).fv(t,i)*seq(n).A(i,j)*seq(n).B(j,seq(n).O(t+1))*seq(n).bv(t+1,j))/(seq(n).fv(t,1)*seq(n).bv(t,1)+seq(n).fv(t,2)*seq(n).bv(t,2)+seq(n).fv(t,3)*seq(n).bv(t,3));
                    seq(n).E(t,i,j)=seq(n).E(t,i,j)/seq(n).Z(t+1);
                    
                end;
                seq(n).R(t,i)=seq(n).E(t,i,1)+seq(n).E(t,i,2)+seq(n).E(t,i,3);
            end;
            
        end;
        for i=1:3
            seq(n).R(l(n),i)=seq(n).fv(l(n),i)/sum(seq(n).fv(l(n),:));
        end;
        
        
        for i=1:3
            sumR=sum(seq(n).R(1:l(n)-1,i));
            for j=1:3
                seq(n).value(i,j)=sum(seq(n).E(:,i,j))/sumR;
                
            end;
        end;
        
        
        
        for i=1:3
            sumR=sum(seq(n).R(:,i));
            
            for j=1:20
                
                seq(n).value2(i,j)=sum(seq(n).R(seq(n).O==j,i))/sumR;
            end;
            
        end;
        
        seq(n).pie2=seq(n).R(1,:);
        seq(n).last=seq(n).current;
        seq(n).current=0;
        for i=1:l(n)
            seq(n).current=seq(n).current+log10(seq(n).Z(i));
        end;
    end;
end;
%calculate state transition probability matrix A, emission probability matrix B and initial state probability distribution PIE
A(1:3,1:3)=0;
B(1:3,1:20)=0;
PIE(1:3)=0;
for n=1:100
    A=A+seq(n).value;
    B=B+seq(n).value2;
    PIE=PIE+seq(n).pie2;
end;
A=A/110;
B=B/110;
PIE=PIE/100;