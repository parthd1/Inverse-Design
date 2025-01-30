function [s11, p1, p2] = modelResult(model)

    pmlBlabels = getOutLabels(model);
    eg1 = model.result.evaluationGroup.create('eg1', 'EvaluationGroup');
    eg1.create('s11', 'EvalGlobal');
    eg1.create('totalPower', 'IntSurface');
    eg1.create('outPower1', 'IntSurface');
    eg1.create('outPower2', 'IntSurface');
    eg1.set('data', 'dset1');
    eg1.set('looplevelinput', {'all'});
    eg1.feature('s11').set('expr', {'abs(solid.Smatrix11)'});
    eg1.feature('s11').set('unit', {'1'});
    eg1.feature('s11').set('descr', {'S11'});
    eg1.feature('s11').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
    %eg1.feature('totalPower').set('expr', {'solid.pzm1.nIpml'});
    %eg1.feature('totalPower').set('unit', {'W'});
    %eg1.feature('totalPower').set('descr', {'Energy flux to surronding PML'});
    %eg1.feature('totalPower').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
    %eg1.feature('totalPower').selection.set([16 21 31 47 pmlBlabels(1)]);
    eg1.feature('outPower1').set('expr', {'solid.pzm1.nIpml'});
    eg1.feature('outPower1').set('unit', {'W'});
    eg1.feature('outPower1').set('descr', {'Energy flux to PML'});
    eg1.feature('outPower1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
    eg1.feature('outPower1').selection.set(pmlBlabels(1));
    eg1.feature('outPower2').set('expr', {'solid.pzm1.nIpml'});
    eg1.feature('outPower2').set('unit', {'W'});
    eg1.feature('outPower2').set('descr', {'Energy flux to PML'});
    eg1.feature('outPower2').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
    eg1.feature('outPower2').selection.set(pmlBlabels(2));
    eg1.run;
    data = mphtable(model, 'eg1');
    s11 = data.data(2);
    p1 = data.data(3);
    p2 = data.data(4);
end

function pmlBlabels = getOutLabels(model)
    para = [model.param.evaluate('L'), model.param.evaluate('outYshift')];
    %coordPmlB1 = [para(1)/2+model.param.evaluate('PML_t');0;model.param.evaluate('AIN_t')/2];
    coordPmlB1 = [para(1)/2+model.param.evaluate('outputL');para(2);model.param.evaluate('AIN_t')/2];
    coordPmlB2 = [para(1)/2+model.param.evaluate('outputL');-para(2);model.param.evaluate('AIN_t')/2];
    pmlBlabels = [mphselectball(model, 'geom1', coordPmlB1, 'boundary', 'radius', 0.01),...
                  mphselectball(model, 'geom1', coordPmlB2, 'boundary', 'radius', 0.01)];
end
