% ================Inference for Individual Predicate based on Thresold:
%clear DATA_OUT_SORT;
%clear SelectTriples;

%t1=thresold1674M;
%t2=thresold1674;
%t3 = t2-t1;
%T = t- t3;
%T = 1;
%AllFactorValue = DATA_OUT(:,5);
Inferred_Entity = DATA_OUT(:,2);
CHECK_DATA_OUT = DATA_OUT;
CHECK_DATA_OUT(:,6) = 0;
%M = Missing_147_60(:,2);
%MissingEntity = GroundTruth_100FilmActor_10324(:,2);

MissingEntity = FilmActor_Missing_List(:,2);
%sum=0;
%n=0;
L = length(Inferred_Entity);
for m=1:L
 tempE = Inferred_Entity(m);
 if (any(MissingEntity(:) == tempE))
     CHECK_DATA_OUT(m, 6)= 1; %change 1 to 2 as for domain relevance is 2.
     
 end
 
end

disp('All rows process completed...');



