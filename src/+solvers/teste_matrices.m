function error = teste_matrices(x,u,idb,family,partialAreas,A,surfaceCorrection,ndof,par_omega,model,damage)
%% 
A = A*ones(length(x));
K_anal = model.analyticalStiffnessMatrix(x,u,ndof,idb,family,partialAreas,surfaceCorrection,A,par_omega,damage,model.history);
K_num = tangentStiffnessMatrix(x,u,idb,family,partialAreas,A,surfaceCorrection,ndof,par_omega,model,damage,model.history);
% Pre-processing
median = prctile(K_anal(K_anal(:)~=0),50);
lower_bound = abs(median)*1e-6;
K_num(abs(K_num) < lower_bound) = 0;
K_anal(abs(K_anal) < lower_bound) = 0;

error = zeros(size(K_anal));
idx = K_anal~= 0;
error(idx) = (K_num(idx) - K_anal(idx))./K_anal(idx);
Comp = [K_num(:) K_anal(:)];
disp(Comp(1:50,:))


end