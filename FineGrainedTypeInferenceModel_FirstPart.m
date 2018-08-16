% Fine-grained Type Inference based on Tensor Factorization & Tensor Reconstruction

clear thresold;
temp=0;
sum=0;
n=0;
p=KG_Type(1,4);
s = KG_Type(1,2);
d =1;
DATA_OUT= [];
%thresold =0;
for o=1:14950
        if (Domain4thOrderTensor(s, o, p, d) == 0)
            temp = Domain4thOrderTensor_OUT(s, o, p, d);
            %temp = temp*(1000000000);
            input01 = [s o p d temp];
            DATA_OUT = [DATA_OUT; input01];
            
            %i = 'Process running: ';
            %X = [i, 'No. ',num2str(n),' select triples: [ ', num2str(s), ' - ', num2str(o), ' - ', num2str(p), ' ] Value: ', num2str(temp)];
            %disp(X);
            
            %disp(o);
           % disp(temp);
            %sum = temp + sum;
            %disp('1000 rows process completed...');
            %disp(s); disp(0);
            n=n+1;
            if (n == 1000)
            disp(n);
            end
            if (n == 2000)
            disp(n);
            end
            if (n == 3000)
            disp(n);
            end
            if (n == 4000)
            disp(n);
            end
            if (n == 5000)
            disp(n);
            end
            if (n == 6000)
            disp(n);
            end
            if (n == 7000)
            disp(n);
            end
            if (n == 9000)
            disp(n);
            end
            if (n == 10000)
            disp(n);
            end
            if (n == 11000)
            disp(n);
            end
              if (n == 12000)
            disp(n);
              end
              if (n == 13000)
            disp(n);
              end
            if (n == 14000)
            disp(n);
            end
            
        end
    %for o=1:100
    
    
    %disp('1000 rows process completed...');
    %n =n+1;
end
%thresold = sum/n;
disp('All rows process completed...');
%disp(sum);
disp(n);
%disp(thresold)
clear p;
clear s;
clear o;
%Thresold_TVActor = thresold;
clear thresold;

% ==============================================

% ================Inference for Individual Predicate based on Thresold:
%clear DATA_OUT_SORT;
%clear SelectTriples;

%t1=thresold1674M;
%t2=thresold1674;
%t3 = t2-t1;
%T = t- t3;
%T = Thresold_FilmActor;
%AllFactorValue = DATA_OUT(:,5);
%M = Missing_147_60(:,2);
%MissingEntity = GroundTruth_100FilmActor_10324(:,2);
%sum=0;
%n=0;
%L = length(AllFactorValue);
%AvgFac =0;
%n=0;
%for m=1:L

 %   FactorValue = AllFactorValue(m);
  %   if (FactorValue > 0)
   %     Sum = Sum + FactorValue;
    %    n = n+1;
     %end
%end

%AvgFac = Sum/n;

%disp(AvgFac);
%Threshold = (T + AvgFac);
%disp(Threshold);
%disp(AvgFac);
%disp(Threshold);
%disp(T);
%nn=0;

%p=32;
%s = 10324;
%d =1;

%LL = length(MissingEntity);
%for mm=1:LL
  %  MissingM = MissingEntity(mm);
  %  tempM = Domain4thOrderTensor_OUT(s, MissingM, p, d);
 % if tempM > Threshold
     %disp('found');
   %  nn = nn+1;
 % end
%end


 %disp('found all');
 %disp(nn);
 %disp('missing all');
 %disp(LL);
 %disp('Threshold');
 %disp(AvgFac);

disp('All rows process completed...');



