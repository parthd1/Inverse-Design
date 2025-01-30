clc;
clear;
addpath('/home/programs/comsol62/mli/');
%addpath('/Applications/COMSOL62/Multiphysics/mli');
import com.comsol.model.*
import com.comsol.model.util.*

try
    mphstart(2036);
   % mphstart('10.149.199.25', 2036, 'local', '31083108!MIM');
catch

end
ModelUtil.remove('Model');
geom_para = [20, 0.4, 0.32, 1, 1]; %geometry parameter n, a, d, w_in, w_out;
coupler_geom = [0.8, 0.5];
param = struct('geom_param', geom_param, 'coupler_param', coupler_param);
fprintf('Geometry parameter: n -- %d; a -- %f um; d -- %f um; w_in -- %f um; w_out -- %f um\n', geom_para(1), geom_para(2), geom_para(3), ...
    geom_para(4), geom_para(5));
splitRatio = [9, 1];
fprintf('Target spliting ratio: %d:%d = %f \n',splitRatio(1), splitRatio(2), splitRatio(1)/splitRatio(2));
hyperparam = 0.2;
fprintf('Chosen hyperparameter: %f\n', hyperparam);
modelUsed = 2; % power spliter model being used: 1 -- no PML boundary; 2 -- PML cover all boundary
if modelUsed == 1
    fprintf('Model being used: 1 (no PML boundary)\n');
elseif modelUsed == 2
    fprintf('Model being used: 2 (PML with full cover)\n');
elseif modelUsed == 3
    fprintf('Model being used: 3 (no PML boundary with coupler)')
elseif modelUsed == 4
    fprintf('Model being used: 4 (PML with full cover with coupler)')
[totalLoss, bestImages] = DBS(param, 0.001, splitRatio, hyperparam, modelUsed);
save('./genData/data0.mat', '-mat');
fprintf('Computation finished!\n');




