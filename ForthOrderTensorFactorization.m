% The 4th Order Tensor Factorization model for fine-grained type inference in knowledge graph and Textual data by ABM Moniruzzaman



% Merging KG and Texual Data
DomainFB215KALL = [DomainFB15K; DomainFB15K_237];
% Constructing 4th Oder Tensor
Domain4thOrderTensor = sptensor(DomainFB215KALL(:,1:4), DomainFB215KALL(:,5));
disp('Construction 4th Oder Tensor Complete'); 
% Factorizing and Reconmstruction of 4th Oder Tensor with cp_opt - CP factorization with Optimization based algorithm. 
Domain4thOrderTensor_OUT = ncp_als(Domain4thOrderTensor, 4);
disp('Factorizing and Reconmstruction of 4th Oder Tensor with cp_opt Complete');
% Here Domain4thOrderTensor_OUT is Reconstructed Tensor 
%=========================================================================