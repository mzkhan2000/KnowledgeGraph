%This function finds Latent Entity Similarity values for Probabilitistic Inference for Fine-graied Type Inference in Knowledge Graph
% Latent Entity Similarity values with 50 rank.


Relations = FB15K_KG_TRAIN_FilmActor(:,3);
SubjectEntity = FB15K_KG_TRAIN_FilmActor(:,1);
ObjectEntity = FB15K_KG_TRAIN_FilmActor(:,2);
%EntryValues = NB_FB15K(:,4);

% populating entity similarity matrix for KG enties for each relations
Matrix_Sub = [SubjectEntity,Relations,ObjectEntity];
Matrix_Obj = [ObjectEntity,Relations,SubjectEntity];
Matrix_KG =  [Matrix_Sub; Matrix_Obj]; 
clear Relations;
clear SubjectEntity;
clear ObjectEntity;
clear Matrix_Sub;
clear Matrix_Obj;

% for delete duplicate entries from entity similarity matrix.
%Matrix_KG = round(Matrix_KG);
%Matrix_KG = unique(Matrix_KG,'rows');

clear NB_FB15K; clear NB_Type; clear Predicate; clear TypedEntity; clear n; clear L; clear tempM; clear m;
% input Knowledge Graph
NB_FB15K = Matrix_KG;
% input type relations
NB_Type = KG_Type(1, 1);
v = KG_Type(1, 3);
% inserting 4th order column in dataset
NB_FB15K(:,4)=1; %change 0 to 1 as 0 is not real positive number.
NB_FB15K(:,5)=1; % Inserting weight to 5th dimention of 4th-order tensor
%Subject = UpdateFB15K(:,1);
%Object = UpdateFB15K(:,2);
% slecting all predicate from KG
Predicate = NB_FB15K(:,2);
%FB15K_INVALID = [Object, Subject, Predicate];
%FilmDomainObject(:,1) = [];
%TypedEntity = [];
n=0;
L = length(Predicate);
for m=1:L
 tempM = Predicate(m);
 if (NB_Type == tempM)
     NB_FB15K(m, 4)= v; % Inserting domain information to 4th dimention of tensor
     NB_FB15K(m, 5)= 2; % Inserting domain information to 4th dimention of tensor
     
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
     n = n+1;
     
 end
 
end
% clear variables used for previous process....
clear NB_Type; clear Predicate; clear TypedEntity; clear n; clear L; clear tempM; clear m; clear v;
clear Matrix_KG;
% 
%Matrix_KG_T = Matrix_KG;
%Matrix_KG_T(:,3) = 1;
%Matrix_KG_T(:,4) = 100000;
%Matrix_KG_T(:,3) = 1;
Matrix_KG_T = sptensor(NB_FB15K(:,1:3), NB_FB15K(:,4));
%ESTO_KG = cp_als(Matrix_KG_T, 1);
%Entity_KG = Fac{1};
[ELS_OUT, G, ELS_KG]   =  ncp_als(Matrix_KG_T, 4);

%[Fac,G,out]   = cmtf_opt(Z,R,'init',init,'alg_options',options); 

% ELS stands for Entity Latent Similarity
ELS = ELS_OUT{1};

%clear Fac; 
%clear Matrix_KG;
%U{1} = Fac{1};
%U{2} = Fac{2};
%U{3} = Fac{3};
%FactorE = U{1};
%U2{2} = Fac{4};
%U2{3} = Fac{5};

%P1 = ktensor(U);
%P2 = ktensor(U2);

%data.RECONSTRUCTED_Tensor01      = P1;
%==========================================================================
% Latent Entity Similarity values for all entities in KG calculate done with 50 rank.
%==========================================================================


%==========================================================================
% Conditional Probabilities Values to Estimate Posterior Probabilities
%==========================================================================
%ELS_OUT(1047)
%------------------finds all entities for this type in KG-----------------
NB_Type = KG_Type(1, 2);
Type_Relation = KG_Type(1, 4);
A= FB15K_KG_TRAIN_FilmActor;
ind1 = A(:,1) == NB_Type;
A1 = A(ind1,:);
ind2 = A1(:,3) == Type_Relation;
A2 = A1(ind2,:);

EntityList = A2(:,2);
EntityList1 = sort(EntityList);
clear A; clear A1; clear A2; clear ind1; clear ind2; clear EntityList; 

n=0;
L = length(EntityList1);
%L = length(ELS);
% addinting extra column index to ELS...
%Index = find(ELS, L);
%ELS = [Index, ELS];
%firstColumn = cell2mat(ELS(:,1));

ELSC1 = [];

for m=1:L
 tempM = EntityList1(m, 1);
 %if (any(EntityList1 == tempM)) % done....
          
    ELS_ROW = ELS(tempM, :);
    ELS_ROW = [tempM, ELS_ROW];
    ELSC1 = [ELSC1; ELS_ROW];
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
     n = n+1;
     
 end
clear n; clear L; clear ELS_ROW; clear tempM;

n=0;
L = length(EntityList1);
%L = length(ELS);
% addinting extra column index to ELS...
%Index = find(ELS, L);
%ELS = [Index, ELS];
%firstColumn = cell2mat(ELS(:,1));
T = KG_Type(1, 1);
ELSC2 = [];

for m=1:L
 tempM = EntityList1(m, 1);
 %if (any(EntityList1 == tempM)) % done....
      ELS_VAL = ELS_OUT(tempM, T);    
     %ELS_ROW = ELS(tempM, :);
    ELSC2 = [ELSC2; ELS_VAL];
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
     n = n+1;
     
 end
clear n; clear L; clear ELS_VAL; clear tempM;

ELSCF = [ELSC1, ELSC2]; 

clear ELSC1; clear ELSC2;

%==========================================================================
% Conditional Probabilities Values to Estimate Posterior Probabilities
% -----------------------NB_NonEntitySimilarity_KG-----------------------
%==========================================================================