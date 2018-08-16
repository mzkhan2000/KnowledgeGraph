
%==========================================================================
% Conditional Probabilities Values to Estimate Posterior Probabilities
%==========================================================================


n=0;

%L = length(ELS);
% addinting extra column index to ELS...
%Index = find(ELS, L);
%ELS = [Index, ELS];
%firstColumn = cell2mat(ELS(:,1));
index = INDEX_ELSC(:,1);
C = setdiff(INDEX_ELSC(:,1),EntityList1);
L = length(C);
ELSC_NON1 = [];
ELSC_NON2 = [];
for m=1:L
    tempM = C(m, 1);
    %if (any(EntityList1(:) == m))
        %any(FilmDomain(:) == tempM)
        ELS_ROW = ELS(tempM, :);
        ELSC_NON1 = [ELSC_NON1; ELS_ROW];
    %end
   
    %tempM = EntityList1(m, 1);
 %if (any(EntityList1 == tempM)) % done....
      ELS_VAL = ELS_OUT(tempM, T);    
     %ELS_ROW = ELS(tempM, :);
    ELSC_NON2 = [ELSC_NON2; ELS_VAL];
    % ELS(1, m);
    %if (any(EntityList1 == tempM)) % done....
          
    
    % ELS(1, m);
    % NB_FB15K(m, 4)= v; % Inserting domain information to 4th dimention of tensor
    % NB_FB15K(m, 5)= 2; % Inserting domain information to 4th dimention of tensor
     
     %Obj = NB_FB15K(m, 2); % Selecting object entities to Sub variable
     %Sub = NB_FB15K(m, 1); % Selecting subject entities to Sub variable 
                % FilmDomainObject(:,1) = Obj;
     %TypedEntity = [TypedEntity; Obj]; %adding object entities vertically to array
     %TypedEntity = [TypedEntity; Sub]; %adding subject entities vertically to array
                %FilmDomainObject = [FilmDomainObject, Obj]; %adding horizontally to array
                %FilmDomainObject = vercat([FilmDomainObject, Obj]);
                %A = cumsum(Obj);
                %FilmDomainObject (:,1) = Obj;
            %  for m=1:L  % this loop start for finding object entities for this domain
       %     tempM = Predicate(m);
        %    if (any(FilmDomain(:) == tempM))
         %       UpdateFB15K(m, 4)= 1;
                % FilmDomainEntity
                % disp('found');
          %      n = n+1;
     
           % end
         %end
    % FilmDomainEntity
    % disp('found');
     m=m+1;
     
 end
clear n; clear L; clear ELS_ROW; clear tempM;

ELSC_NON = [ELSC_NON1, ELSC_NON2];

clear ELSC_NON1; clear ELSC_NON2;
