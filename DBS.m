function [totalLoss, bestImages] = DBS(geom_params, sigma, target, hyperparam, modelUsed)
    rng("shuffle");
    %binaryImage = randi(2, geom_param(1))>1;
    binaryImage = get_init_image('init_best.mat');
    fprintf('The total holes for random initilization is : %d\n', sum(binaryImage, 'all'))
    iterNum = 1;
    totalLoss = [];
    bestImages = [];
    [loss, bestImage] = random_iter(binaryImage, geom_params, target, hyperparam, modelUsed, iterNum);
    totalLoss = [totalLoss, loss];
    bestImages = [bestImages, bestImage];
    
    while loss > sigma
        iterNum = iterNum + 1;
        [loss, bestImage] = random_iter(bestImage, geom_params, target, hyperparam, modelUsed, iterNum);
        totalLoss = [totalLoss, loss];
        bestImages = cat(3, bestImages, bestImage);
    end
    
    function init_image = get_init_image(name)
        init_image = load(name).best;
    end

   % function rm_plot = generate_rm_plot(n)
   %     rm_plot = randi(2, n);
   %     rm_plot = rm_plot > 1;
   % end
end





