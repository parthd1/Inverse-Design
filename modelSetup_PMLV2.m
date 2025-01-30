%%% this modelSetup function create model using boundary load as input 
function output = modelSetup_PMLV2(model, binaryImage, geom_params)
%%%%%%%%%%%%%%%%%%%%%%%%%% Set up model parameter %%%%%%%%%%%%%%%%%%%%%%%%%
    geom_param = geom_params.geom_param;
    n = geom_param(1); % number of cells
    a = geom_param(2); % um;0.23
    d = geom_param(3); % um;radius of hole;0.18
    L = a *n;
    %w = 1.5; % um;waveguide width
    w_in = geom_param(4);
    w_out = geom_param(5);
    %outYshift = L/2 - w/2; % output unit y-axis shift
    outYshift = L/2 - w_out/2;
    model.param.set('a', string(a));
    model.param.descr('a', 'cell length');
    model.param.set('n', string(n));
    model.param.descr('n', 'number of pixiel');
    model.param.set('L', string(L));
    model.param.descr('L', 'Length of the scattering pattern');
    model.param.set('inputL', '5');
    model.param.descr('inputL', 'Length of the input section');
    model.param.set('outputL', '3');
    model.param.descr('outputL', 'ength of the output section');
    model.param.set('outYshift', string(outYshift));
    model.param.descr('outYshift', 'y-axis shift for output waveguide');
    %model.param.set('w', string(w));
    %model.param.descr('w', 'width of In/Output waveguide');
    model.param.set('w_out', string(w_out));
    model.param.descr('w_out', 'width of output waveguide');
    model.param.set('w_in', string(w_in));
    model.param.descr('w_in', 'width of input waveguide');
    model.param.set('x_out1', 'L/2 + 2');
    model.param.set('x_out2', 'x_out1');
    model.param.set('y_out1', 'L/4');
    model.param.set('y_out2', '-L/4');
    model.param.descr('x_out1', 'x-cord of the output section 1');
    model.param.descr('y_out1', 'y-cord of the output section 1');
    model.param.descr('x_out2', 'x-cord of the output section 2');
    model.param.descr('y_out2', 'y-cord of the otput section 3');
    model.param.set('x_c', '0');
    model.param.set('y_c', '0');
    model.param.descr('x_c', 'x-cord of center');
    model.param.descr('y_c', 'y-cord of center');
    model.param.set('d', string(d));
    model.param.set('AIN_t', '0.2');
    model.param.descr('d', 'diameter of air hole');
    model.param.descr('AIN_t', 'thisckness of LN');
    model.param.set('ka', '3.37*1e6');
    model.param.descr('ka', 'A0 mode wave number');
    model.param.set('wl', '2 * pi/ka * 1e6');
    model.param.descr('wl', 'A0 mode wave length');
    model.param.set('PML_t', 'wl');
    model.param.descr('PML_t', 'thickness of PML');
    model.param.set('etchSep', '0.1');
    model.param.descr('etchSep', 'seperation dist for etching');
    model.param.set('buffer_t', 'wl/2');
    model.param.descr('buufer_t', 'input PML buffer layer thickness');
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Create model component %%%%%%%%%%%%%%%%%%%%%%%%%%

    comp1 = model.component.create('comp1', true);
    
    geom1 = comp1.geom.create('geom1', 3);
    geom1.lengthUnit([native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']); % set unit for geom1
    mesh1 = comp1.mesh.create('mesh1');
    
    geom1.create('wp1', 'WorkPlane');
    
    geom1.feature('wp1').geom.create('pol1', 'Polygon');
    geom1.feature('wp1').geom.feature('pol1').set('source', 'table');
    geom1.feature('wp1').geom.feature('pol1').set('table', {'x_c + -L/2' 'y_c + L/2';  ...
        'x_c - L/2' 'y_c + w_in/2';  ...
        'x_c - L/2 - inputL' 'y_c + w_in/2';  ...
        'x_c - L/2 - inputL' 'y_c - w_in/2';  ...
        'x_c - L/2' 'y_c - w_in/2';  ...
        'x_c - L/2' 'y_c - L/2';  ...
        'x_c + L/2' 'y_c - L/2';  ...
        'x_c + L/2' 'y_c - outYshift - w_out/2';  ...
        'x_c + L/2 + outputL' 'y_c - outYshift - w_out/2';  ...
        'x_c + L/2 + outputL' 'y_c - outYshift + w_out/2';  ...
        'x_c + L/2' 'y_c - outYshift + w_out/2';  ...
        'x_c + L/2' 'y_c + outYshift - w_out/2';  ...
        'x_c + L/2 + outputL' 'y_c + outYshift - w_out/2'; ...
        'x_c + L/2 + outputL' 'y_c + outYshift + w_out/2'; ...
        'x_c + L/2' 'y_c + outYshift + w_out/2'; ...
        'x_c + L/2' 'y_c + L/2';});
    geom1.feature('wp1').geom.feature('pol1').set('type', 'solid');
    
    geom1.feature('wp1').geom.create('c1', 'Circle');
    geom1.feature('wp1').geom.feature('c1').set('r', 'd/2');
    geom1.feature('wp1').geom.feature('c1').set('pos', {'-L/2 + a/2','L/2 - a/2'});

    geom1.feature('wp1').geom.create('copy1', 'Copy');
    geom1.feature('wp1').geom.feature('copy1').selection('input').init;
    geom1.feature('wp1').geom.feature('copy1').selection('input').set({'c1'});

    disp = plot2disp(binaryImage, a);
    dispx = disp2str(disp(1, :));
    dispy = disp2str(disp(2, :));
    geom1.feature('wp1').geom.feature('copy1').set('displx',dispx);
    geom1.feature('wp1').geom.feature('copy1').set('disply',dispy);
    geom1.feature('wp1').geom.feature('copy1').set('keep', false);
    geom1.feature('wp1').geom.feature.create('diff1', 'Difference');
    geom1.feature('wp1').geom.feature('diff1').selection('input').set({'pol1'});
    geom1.feature('wp1').geom.feature('diff1').selection('input2').set({'copy1'});
    geom1.feature('wp1').label('Scattering part 2D');
    
    geom1.feature('wp1').geom.create('ls1', 'LineSegment');
    geom1.feature('wp1').geom.feature('ls1').set('specify1', 'coord');
    geom1.feature('wp1').geom.feature('ls1').set('coord1', {'-L/2' 'w_in/2'});
    geom1.feature('wp1').geom.feature('ls1').set('specify2', 'coord');
    geom1.feature('wp1').geom.feature('ls1').set('coord2', {'-L/2' '-w_in/2'});
    geom1.feature('wp1').geom.create('ls2', 'LineSegment');
    geom1.feature('wp1').geom.feature('ls2').set('specify1', 'coord');
    geom1.feature('wp1').geom.feature('ls2').set('coord1', {'L/2' 'w_out/2 + outYshift'});
    geom1.feature('wp1').geom.feature('ls2').set('specify2', 'coord');
    geom1.feature('wp1').geom.feature('ls2').set('coord2', {'L/2' '-w_out/2 + outYshift'});
    geom1.feature('wp1').geom.create('ls3', 'LineSegment');
    geom1.feature('wp1').geom.feature('ls3').set('specify1', 'coord');
    geom1.feature('wp1').geom.feature('ls3').set('coord1', {'L/2' 'w_out/2 - outYshift'});
    geom1.feature('wp1').geom.feature('ls3').set('specify2', 'coord');
    geom1.feature('wp1').geom.feature('ls3').set('coord2', {'L/2' '-w_out/2 - outYshift'});
    
    geom1.feature('wp1').geom.create('outer region1', 'Rectangle');
    geom1.feature('wp1').geom.feature('outer region1').set('size', {'L' 'PML_t'});
    geom1.feature('wp1').geom.feature('outer region1').set('pos', {'-L/2' 'L/2'});
    geom1.feature('wp1').geom.create('outer region2', 'Rectangle');
    geom1.feature('wp1').geom.feature('outer region2').set('size', {'L' 'PML_t'});
    geom1.feature('wp1').geom.feature('outer region2').set('pos', {'-L/2' '-L/2 - PML_t'});
    geom1.feature('wp1').geom.create('outer region3', 'Rectangle');
    geom1.feature('wp1').geom.feature('outer region3').set('size', {'PML_t' 'L/2 - w_in/2 - etchSep'});
    geom1.feature('wp1').geom.feature('outer region3').set('pos', {'-L/2 - PML_t' 'w_in/2 + etchSep'});
    geom1.feature('wp1').geom.create('outer region4', 'Rectangle');
    geom1.feature('wp1').geom.feature('outer region4').set('size', {'PML_t' 'L/2 - w_in/2 - etchSep'});
    geom1.feature('wp1').geom.feature('outer region4').set('pos', {'-L/2 - PML_t' '-L/2'});
    geom1.feature('wp1').geom.create('outer region5', 'Rectangle');
    geom1.feature('wp1').geom.feature('outer region5').set('size', {'PML_t' 'L - w_out*2 - etchSep*2'});
    geom1.feature('wp1').geom.feature('outer region5').set('pos', {'L/2' '-L/2 + w_out + etchSep'});
    
    geom1.feature('wp1').geom.create('PML1', 'Rectangle');
    geom1.feature('wp1').geom.feature('PML1').set('size', {'L' 'PML_t'});
    geom1.feature('wp1').geom.feature('PML1').set('pos', {'-L/2' 'L/2 + PML_t'});
    geom1.feature('wp1').geom.create('PML2', 'Rectangle');
    geom1.feature('wp1').geom.feature('PML2').set('size', {'L' 'PML_t'});
    geom1.feature('wp1').geom.feature('PML2').set('pos', {'-L/2' '-L/2 - PML_t - PML_t'});
    geom1.feature('wp1').geom.create('PML3', 'Rectangle');
    geom1.feature('wp1').geom.feature('PML3').set('size', {'PML_t' 'L/2 - w_in/2 - etchSep'});
    geom1.feature('wp1').geom.feature('PML3').set('pos', {'-L/2 - PML_t*2' 'w_in/2 + etchSep'});
    geom1.feature('wp1').geom.create('PML4', 'Rectangle');
    geom1.feature('wp1').geom.feature('PML4').set('size', {'PML_t' 'L/2 - w_in/2 - etchSep'});
    geom1.feature('wp1').geom.feature('PML4').set('pos', {'-L/2 - PML_t*2' '-L/2'});
    geom1.feature('wp1').geom.create('PML5', 'Rectangle');
    geom1.feature('wp1').geom.feature('PML5').set('size', {'PML_t' 'L - w_out*2 - etchSep*2'});
    geom1.feature('wp1').geom.feature('PML5').set('pos', {'L/2 + PML_t' '-L/2 + w_out + etchSep'});
    geom1.feature('wp1').geom.create('outPML1', 'Rectangle');
    geom1.feature('wp1').geom.feature('outPML1').set('size', {'PML_t' 'w_out'});
    geom1.feature('wp1').geom.feature('outPML1').set('pos', {'outputL + L/2' 'outYshift - w_out/2'});
    geom1.feature('wp1').geom.create('outPML2', 'Rectangle');
    geom1.feature('wp1').geom.feature('outPML2').set('size', {'PML_t' 'w_out'});
    geom1.feature('wp1').geom.feature('outPML2').set('pos', {'outputL + L/2' '-outYshift - w_out/2'});
    % Create input geom 
    geom1.feature.create('ext1', 'Extrude');
    geom1.feature('ext1').selection('input').set({'wp1'});
    geom1.feature('ext1').setIndex('distance', 'AIN_t', 0);
    geom1.run;
    %%%%%%%%%%%%%%%%%%%%%%%%% Set up PML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    comp1.coordSystem.create('pml1', 'PML');
    comp1.coordSystem('pml1').selection.set([2 3 6 10 14 15 16]);
    
    %%%%%%%%%%%%%%%%%%%%%%% Create material %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% AIN %%%%%%%%%%%%%%%%%
    AIN = comp1.material.create('mat1', 'Common');
    AIN.propertyGroup.create('StrainCharge', 'Strain-charge form');
    AIN.propertyGroup.create('StressCharge', 'Stress-charge form');
    AIN.label('Aluminum Nitride');
    AIN.set('family', 'custom');
    AIN.set('customspecular', [0.7843137254901961 1 1]);
    AIN.set('diffuse', 'custom');
    AIN.set('customdiffuse', [0.7843137254901961 0.7843137254901961 0.7843137254901961]);
    AIN.set('ambient', 'custom');
    AIN.set('customambient', [0.7843137254901961 0.7843137254901961 0.7843137254901961]);
    AIN.set('noise', true);
    AIN.set('fresnel', 0.9);
    AIN.set('roughness', 0.1);
    AIN.set('diffusewrap', 0);
    AIN.set('reflectance', 0);
    AIN.propertyGroup('def').set('relpermittivity', {'9' '0' '0' '0' '9' '0' '0' '0' '9'});
    AIN.propertyGroup('def').set('density', '3300[kg/m^3]');
    AIN.propertyGroup('StrainCharge').set('sE', {'2.8987e-12[1/Pa]' '-9.326e-013[1/Pa]' '-5.0038e-013[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '-9.326e-013[1/Pa]' '2.8987e-12[1/Pa]' '-5.0038e-013[1/Pa]' '0[1/Pa]'  ...
        '0[1/Pa]' '0[1/Pa]' '-5.0038e-013[1/Pa]' '-5.0038e-013[1/Pa]' '2.8253e-012[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]'  ...
        '0[1/Pa]' '8E-12[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '8E-12[1/Pa]' '0[1/Pa]'  ...
        '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '0[1/Pa]' '7.6628E-12[1/Pa]'});
    AIN.propertyGroup('StrainCharge').set('dET', {'0[C/N]' '0[C/N]' '-1.9159e-012[C/N]' '0[C/N]' '0[C/N]' '-1.9159e-012[C/N]' '0[C/N]' '0[C/N]' '4.9597e-012[C/N]' '0[C/N]'  ...
        '-3.84e-012[C/N]' '0[C/N]' '-3.84e-012[C/N]' '0[C/N]' '0[C/N]' '0[C/N]' '0[C/N]' '0[C/N]'});
    AIN.propertyGroup('StrainCharge').set('epsilonrT', {'9.2081' '0' '0' '0' '9.2081' '0' '0' '0' '10.1192'});
    AIN.propertyGroup('StressCharge').set('cE', {'4.1e+011[Pa]' '1.49e+011[Pa]' '0.99e+011[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '1.49e+011[Pa]' '4.1e+011[Pa]' '0.99e+011[Pa]' '0[Pa]'  ...
        '0[Pa]' '0[Pa]' '0.99e+011[Pa]' '0.99e+011[Pa]' '3.89e+011[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '0[Pa]'  ...
        '0[Pa]' '1.25e+011[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '1.25e+011[Pa]' '0[Pa]'  ...
        '0[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '0[Pa]' '1.305e+011[Pa]'});
    AIN.propertyGroup('StressCharge').set('eES', {'0[C/m^2]' '0[C/m^2]' '-0.58[C/m^2]' '0[C/m^2]' '0[C/m^2]' '-0.58[C/m^2]' '0[C/m^2]' '0[C/m^2]' '1.55[C/m^2]' '0[C/m^2]'  ...
        '-0.48[C/m^2]' '0[C/m^2]' '-0.48[C/m^2]' '0[C/m^2]' '0[C/m^2]' '0[C/m^2]' '0[C/m^2]' '0[C/m^2]'});
    AIN.propertyGroup('StressCharge').set('epsilonrS', {'9' '0' '0' '0' '9' '0' '0' '0' '9'});
    
    %%%%%%%%%%%%%%%%%%%%%% Set up Physics %%%%%%%%%%%%%%%%%%%%%%
    comp1.physics.create('solid', 'SolidMechanics', 'geom1');
    comp1.physics('solid').create('pzm1', 'PiezoelectricMaterialModel', 3);
    comp1.physics('solid').feature('pzm1').selection.set([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);
    comp1.physics('solid').create('port1', 'Port', 2);
    comp1.physics('solid').feature('port1').selection.set([1]);
    comp1.physics('solid').feature('port1').set('power', 0.01);
    
    %%%%%%%%%%%%%%%%%%% Set up mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [outElabels, pmlElabels, outBlabels, pmlMeshBlabels] = getMeshLables(model);
    mesh1.create('ftri1', 'FreeTri');%mesh for scattering
    mesh1.feature('ftri1').selection.geom('geom1', 2);
    mesh1.feature('ftri1').selection.set([37]);
    mesh1.feature('ftri1').create('size1', 'Size');
    mesh1.feature('ftri1').feature('size1').set('hauto', 3);

    mesh1.create('map1', 'Map');%mesh for waveguide
    mesh1.feature('map1').selection.geom('geom1', 2);
    mesh1.feature('map1').selection.set([4 outBlabels(1) outBlabels(2)]);
    mesh1.feature('map1').create('dis1', 'Distribution');
    mesh1.feature('map1').feature('dis1').selection.set([8]);
    mesh1.feature('map1').feature('dis1').set('numelem', 20);
    mesh1.feature('map1').create('dis2', 'Distribution');
    mesh1.feature('map1').feature('dis2').selection.set([outElabels(3) outElabels(4)]);
    mesh1.feature('map1').feature('dis2').set('numelem', 12);

    mesh1.create('map2', 'Map');
    mesh1.feature('map2').selection.geom('geom1', 2);
    mesh1.feature('map2').selection.set([19 24 33 45 pmlMeshBlabels(1)]);
    mesh1.feature('map2').create('dis1', 'Distribution');
    mesh1.feature('map2').feature('dis1').selection.set([29 40 49 71 pmlElabels(4)]);
    mesh1.feature('map2').feature('dis1').set('numelem', 4);

    mesh1.create('bl1', 'BndLayer');%mesh for waveguide boundary
    mesh1.feature('bl1').selection.geom('geom1', 2);
    mesh1.feature('bl1').selection.set([4 outBlabels(1) outBlabels(2)]);
    mesh1.feature('bl1').create('blp', 'BndLayerProp');
    mesh1.feature('bl1').feature('blp').selection.set([61 outElabels(1) outElabels(2)]);
    mesh1.feature('bl1').feature('blp').set('blnlayers', 4);
    mesh1.feature('bl1').feature('blp').set('blstretch', 1.5);
    mesh1.feature('bl1').feature('blp').set('blhminfact', 1.5);

    mesh1.create('bl2', 'BndLayer');%mesh for scattering PML boundary
    mesh1.feature('bl2').selection.geom('geom1', 2);
    mesh1.feature('bl2').selection.set([9 14 29 49 pmlMeshBlabels(2)]);
    mesh1.feature('bl2').create('blp', 'BndLayerProp');
    mesh1.feature('bl2').feature('blp').selection.set([28 36 50 77 pmlElabels(1)]);
    mesh1.feature('bl2').feature('blp').set('blnlayers', 3);
    mesh1.feature('bl2').feature('blp').set('blstretch', 1.5);
    mesh1.feature('bl2').feature('blp').set('blhminfact', 2);

    mesh1.create('bl3', 'BndLayer');%mesh for wavguide PML boundary
    mesh1.feature('bl3').selection.geom('geom1', 2);
    mesh1.feature('bl3').selection.set(pmlMeshBlabels(3:end));
    mesh1.feature('bl3').create('blp', 'BndLayerProp');
    mesh1.feature('bl3').feature('blp').selection.set(pmlElabels(2:end));
    mesh1.feature('bl3').feature('blp').set('blnlayers', 3);
    mesh1.feature('bl3').feature('blp').set('blstretch', 1.5);
    mesh1.feature('bl3').feature('blp').set('blhminfact', 1.5);

    mesh1.create('swe1', 'Sweep');% sweep the mesh
    mesh1.feature('swe1').selection.geom('geom1', 3);
    mesh1.feature('swe1').selection.add([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);
    mesh1.feature('swe1').create('dis1', 'Distribution');
    mesh1.feature('swe1').feature('dis1').set('numelem', 3);
    
    mesh1.run
    output  = model;
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Set up Stuady %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

%mphsave(model, 'powerDemulti_PML_test.mph')



function [outElabels, pmlElabels, outBlabels, pmlMeshBlabels] = getMeshLables(model)
    para = [model.param.evaluate('L'), model.param.evaluate('outYshift'), model.param.evaluate('w_out')];
    coordout1 = [para(1)/2;para(2);model.param.evaluate('AIN_t')];
    coordout2 = [para(1)/2;-para(2);model.param.evaluate('AIN_t')];
    coordoutU1 = [para(1)/2+1;para(2)+para(3)/2;model.param.evaluate('AIN_t')];
    coordoutU2 = [para(1)/2+1;-para(2)+para(3)/2;model.param.evaluate('AIN_t')];
    outElabels = [mphselectball(model, 'geom1', coordout1, 'edge', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordout2, 'edge', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordoutU1, 'edge', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordoutU2, 'edge', 'radius', 0.01)];
    coordPML1 = [para(1)/2+model.param.evaluate('PML_t');0;model.param.evaluate('AIN_t')];
    coordPML2 = [para(1)/2+model.param.evaluate('outputL');para(2);model.param.evaluate('AIN_t')];
    coordPML3 = [para(1)/2+model.param.evaluate('outputL');-para(2);model.param.evaluate('AIN_t')];
    coordPML4 = [para(1)/2+model.param.evaluate('PML_t')/2;para(1)/2-para(3)-model.param.evaluate('etchSep');model.param.evaluate('AIN_t')];
    pmlElabels = [mphselectball(model, 'geom1', coordPML1, 'edge', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordPML2, 'edge', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordPML3, 'edge', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordPML4, 'edge', 'radius', 0.01)];
    coordoutB1 = [para(1)/2 + 1;para(2);model.param.evaluate('AIN_t')];
    coordoutB2 = [para(1)/2 + 1;-para(2);model.param.evaluate('AIN_t')];
    outBlabels = [mphselectball(model, 'geom1', coordoutB1, 'boundary', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordoutB2, 'boundary', 'radius', 0.01)];
    coordPmlMeshB1 = [para(1)/2 + model.param.evaluate('PML_t')/2;0;model.param.evaluate('AIN_t')];
    coordPmlMeshB2 = [para(1)/2 + model.param.evaluate('PML_t')*3/2;0;model.param.evaluate('AIN_t')];
    coordPmlMeshB3 = [para(1)/2+model.param.evaluate('outputL')+model.param.evaluate('PML_t')/2;para(2);model.param.evaluate('AIN_t')];
    coordPmlMeshB4 = [para(1)/2+model.param.evaluate('outputL')+model.param.evaluate('PML_t')/2;-para(2);model.param.evaluate('AIN_t')];
    pmlMeshBlabels = [mphselectball(model, 'geom1', coordPmlMeshB1, 'boundary', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordPmlMeshB2, 'boundary', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordPmlMeshB3, 'boundary', 'radius', 0.01),...
        mphselectball(model, 'geom1', coordPmlMeshB4, 'boundary', 'radius', 0.01)];

end

function disp = plot2disp(plot, a)
    n = length(plot);
    xcord = ones(n) .* [0 : n-1];
    ycord = -ones(n) .* [0 : n-1]';
    dispx = plot .* xcord * a;
    dispy = plot .* ycord * a;
    dispx(find(plot == 0)) = [];
    dispy(find(plot == 0)) = [];
    disp = [dispx;dispy];
end

function disp_str= disp2str(disp)
    str = string(disp);
    disp_str = '';
    for i = 1:length(disp) - 1
        disp_str  = append(disp_str, str(i), ', ');
    end
    disp_str = append(disp_str, str(end));
end


