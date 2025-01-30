function [loss, p1, p2, s11] = update(model, binaryImage, geom_params, target, hyperparam, modelUsed)
    %%%%%%%%%%%%%%%%%%%%%%%%%% Set up model parameter %%%%%%%%%%%%%%%%%%%%%%%%%
    if modelUsed == 1
        model = modelSetup_noPML(model, binaryImage, geom_params);
    elseif modelUsed == 2
        model = modelSetup_PML(model, binaryImage, geom_params);
    elseif modelUsed == 3
        model = modelSetup_PML_coupler(model, binaryImage, geom_params);
    elseif modelUsed == 4
        model = modelSteup_noPML_coupler(model, binaryImage, geom_params);
    model = runSol(model);
    [s11, p1, p2] = modelResult(model);
    gamma1 = target(1)/(target(1) + target(2));
    gamma2 = target(2)/(target(1) + target(2));
    loss = (1 - hyperparam) * (abs(p1 - gamma1 * 0.01) + abs(p2 - gamma2 * 0.01))/0.01 + hyperparam*(s11);
end
