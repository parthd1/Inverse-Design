function [totalLoss, bestImage] = random_iter(binaryImage, geom_params, target, hyperparam, modelUsed, iterNum)
    import com.comsol.model.*
    import com.comsol.model.util.*
    fprintf('Current device iteration: %d \n', iterNum);
    % initilizaing storaged data
    n = geom_params.geom_param(1);
    binaryImages = ones(n, n, n*n+1);
    power1 = ones(1, n*n+1);
    power2 = ones(1, n*n+1);
    sMatrix = ones(1, n*n+1);
    losses = ones(1, n*n+1);
    iter = 1;
    % initilizaing model
    model = ModelUtil.create('Model');
    [lossi, p1, p2, s11] = update(model, binaryImage, geom_params, target, hyperparam, modelUsed);
    binaryImages(:, :, 1) = binaryImage;
    bestImage = binaryImage;
    fprintf('Current iteration: %d -- Loss: %f -- Ratio: %f -- S11: %f -- P1: %f -- P2: %f\n', iter, lossi, p1/p2, s11, p1, p2);
    power1(1, 1) = p1;
    power2(1, 1) = p2;
    sMatrix(1, 1) = s11;
    losses(1, 1) = lossi;
    % iterate randomly all the pixiels
    colIndex = repmat(1:n, 1, n);
    rowIndex = repmat(1:n, 1, n);
    lossf = lossi;
    while size(colIndex, 2) >= 1
        iter = iter + 1;
        model.sol('sol1').clearSolution;
        ModelUtil.remove('Model');
        model = ModelUtil.create('Model');
        % update image and col/row index
        [binaryImage, rowIndex, colIndex] = next_iter(bestImage, rowIndex, colIndex);
        % update simulation
        [loss, p1, p2, s11] = update(model, binaryImage, geom_params, target, hyperparam, modelUsed);
        fprintf('Current iteration: %d -- Loss: %f -- Ratio: %f -- S11: %f -- P1: %f -- P2: %f\n', iter, loss, p1/p2, s11, p1, p2);
        binaryImages(:, :, iter) = binaryImage;
        power1(1, iter) = p1;
        power2(1, iter) = p2;
        sMatrix(1, iter) = s11; % binaryImages, powers, and sMatrix are saved for all calcualted simulstions
        % compared and update lost
        if loss < lossf
            lossf = loss;
            bestImage = binaryImage;
        end
        losses(1, iter) = lossf; % loss function value are saved only for updated loss value
    end
    totalLoss = lossi - lossf;
    save("./genData/data"+int2str(iterNum)+".mat", 'binaryImages', 'bestImage', 'power1', 'power2', 'sMatrix', 'losses');
    model.sol('sol1').clearSolution;
    ModelUtil.remove('Model');
    fprintf('Total improvemnt: %f \n', totalLoss);

       


    function [imgNext, rowIndex, colIndex] = next_iter(binaryImage, rowIndex, colIndex)
        l = size(rowIndex, 2);
        ri = randi(l);
        ci = randi(l);
        r = rowIndex(ri);
        c = colIndex(ci);
        rowIndex(ri) = [];
        colIndex(ci) = [];
        if binaryImage(r, c) == 1
            binaryImage(r, c) = 0;
        else
            binaryImage(r, c) = 1;
        end
        imgNext = binaryImage;

    end


end

