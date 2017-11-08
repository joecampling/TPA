function out = TPA3(n_Si, n_SiO2, lam)
%
% TPA3.m
%
% Model exported on Nov 8 2017, 11:17 by COMSOL 5.3.0.260.

import com.comsol.model.*
import com.comsol.model.util.*

TPA3 = ModelUtil.create('Model');

TPA3.modelPath('C:\Users\jc3e16\Documents\TPA a-Si-H waveguide');

TPA3.label('a-Si-H waveguide test 2_2.mph');

TPA3.comments(['Untitled\n\n']);

TPA3.component.create('comp1', true);

TPA3.component('comp1').geom.create('geom1', 2);

TPA3.result.table.create('evl2', 'Table');
TPA3.result.table.create('tbl1', 'Table');

TPA3.component('comp1').mesh.create('mesh1');

TPA3.component('comp1').geom('geom1').lengthUnit('nm');
TPA3.component('comp1').geom('geom1').scaleUnitValue(true);
TPA3.component('comp1').geom('geom1').create('r1', 'Rectangle');
TPA3.component('comp1').geom('geom1').feature('r1').label('Waveguide');
TPA3.component('comp1').geom('geom1').feature('r1').set('pos', [-250 0]);
TPA3.component('comp1').geom('geom1').feature('r1').set('size', [500 220]);
TPA3.component('comp1').geom('geom1').create('r2', 'Rectangle');
TPA3.component('comp1').geom('geom1').feature('r2').label('Substrate');
TPA3.component('comp1').geom('geom1').feature('r2').set('pos', [-5000 -1950]);
TPA3.component('comp1').geom('geom1').feature('r2').set('size', {'10000' '1950.0'});
TPA3.component('comp1').geom('geom1').create('r7', 'Rectangle');
TPA3.component('comp1').geom('geom1').feature('r7').label('PML (substrate)');
TPA3.component('comp1').geom('geom1').feature('r7').set('pos', [-5500 -2450]);
TPA3.component('comp1').geom('geom1').feature('r7').set('size', [11000 2450]);
TPA3.component('comp1').geom('geom1').create('r3', 'Rectangle');
TPA3.component('comp1').geom('geom1').feature('r3').label('Air');
TPA3.component('comp1').geom('geom1').feature('r3').set('pos', [-5000 0]);
TPA3.component('comp1').geom('geom1').feature('r3').set('size', {'10000.0' '5000.0'});
TPA3.component('comp1').geom('geom1').create('r4', 'Rectangle');
TPA3.component('comp1').geom('geom1').feature('r4').label('PML (air)');
TPA3.component('comp1').geom('geom1').feature('r4').set('pos', [-5500 0]);
TPA3.component('comp1').geom('geom1').feature('r4').set('size', {'11000.0' '5500'});
TPA3.component('comp1').geom('geom1').create('r5', 'Rectangle');
TPA3.component('comp1').geom('geom1').feature('r5').active(false);
TPA3.component('comp1').geom('geom1').feature('r5').label('Silicon base');
TPA3.component('comp1').geom('geom1').feature('r5').set('pos', [-5000 -4000]);
TPA3.component('comp1').geom('geom1').feature('r5').set('size', [10000 2050]);
TPA3.component('comp1').geom('geom1').create('r6', 'Rectangle');
TPA3.component('comp1').geom('geom1').feature('r6').active(false);
TPA3.component('comp1').geom('geom1').feature('r6').label('PML (silicon base)');
TPA3.component('comp1').geom('geom1').feature('r6').set('pos', [-5500 -4500]);
TPA3.component('comp1').geom('geom1').feature('r6').set('size', [11000 2550]);
TPA3.component('comp1').geom('geom1').run;

TPA3.component('comp1').material.create('mat1', 'Common');
TPA3.component('comp1').material.create('mat2', 'Common');
TPA3.component('comp1').material.create('mat3', 'Common');
TPA3.component('comp1').material('mat1').selection.set([5]);
TPA3.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
TPA3.component('comp1').material('mat2').selection.set([1 3]);
TPA3.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
TPA3.component('comp1').material('mat3').selection.set([2 4]);
TPA3.component('comp1').material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');

TPA3.component('comp1').coordSystem.create('pml1', 'PML');
TPA3.component('comp1').coordSystem('pml1').selection.set([1 2]);

TPA3.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
TPA3.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 1);
TPA3.component('comp1').physics('ewfd').feature('sctr1').selection.set([8 10 17]);

TPA3.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
TPA3.component('comp1').mesh('mesh1').create('ref1', 'Refine');

TPA3.result.table('evl2').label('Evaluation 2D');
TPA3.result.table('evl2').comments('Interactive 2D values');
TPA3.result.table('tbl1').comments('Surface Integration 1 (ewfd.normE^2, sqrt(realdot(ewfd.Ex,ewfd.Ex)+realdot(ewfd.Ey,ewfd.Ey))^2)');

TPA3.component('comp1').view('view1').axis.set('xmin', -6207.115234375);
TPA3.component('comp1').view('view1').axis.set('xmax', 6207.115234375);
TPA3.component('comp1').view('view1').axis.set('ymin', -2847.5);
TPA3.component('comp1').view('view1').axis.set('ymax', 5897.5);
TPA3.component('comp1').view('view1').axis.set('abstractviewlratio', -0.06428319960832596);
TPA3.component('comp1').view('view1').axis.set('abstractviewrratio', 0.06428319960832596);
TPA3.component('comp1').view('view1').axis.set('abstractviewbratio', -0.05000000074505806);
TPA3.component('comp1').view('view1').axis.set('abstractviewtratio', 0.05000000074505806);
TPA3.component('comp1').view('view1').axis.set('abstractviewxscale', 15.288460731506348);
TPA3.component('comp1').view('view1').axis.set('abstractviewyscale', 15.288461685180664);

TPA3.component('comp1').material('mat1').label('Silicon');
TPA3.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
TPA3.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
TPA3.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', n_Si);
TPA3.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
TPA3.component('comp1').material('mat2').label('Silica');
TPA3.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', '');
TPA3.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('ki', '');
TPA3.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', n_SiO2);
TPA3.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
TPA3.component('comp1').material('mat3').label('Air');
TPA3.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', '');
TPA3.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', '');
TPA3.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
TPA3.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});

TPA3.component('comp1').coordSystem('pml1').set('ScalingType', 'userDefined');
TPA3.component('comp1').coordSystem('pml1').set('directions', '2');
TPA3.component('comp1').coordSystem('pml1').set('dmax', {'1e-12' '1e-12' '1'});

TPA3.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
TPA3.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
TPA3.component('comp1').mesh('mesh1').feature('size').set('hmax', 100);
TPA3.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.279);
TPA3.component('comp1').mesh('mesh1').feature('ref1').set('boxcoord', true);
TPA3.component('comp1').mesh('mesh1').feature('ref1').set('xmax', 1000);
TPA3.component('comp1').mesh('mesh1').feature('ref1').set('xmin', -1000);
TPA3.component('comp1').mesh('mesh1').feature('ref1').set('ymax', 500);
TPA3.component('comp1').mesh('mesh1').feature('ref1').set('ymin', -1000);
TPA3.component('comp1').mesh('mesh1').run;

TPA3.study.create('std1');
TPA3.study('std1').create('mode', 'ModeAnalysis');

TPA3.sol.create('sol1');
TPA3.sol('sol1').study('std1');
TPA3.sol('sol1').attach('std1');
TPA3.sol('sol1').create('st1', 'StudyStep');
TPA3.sol('sol1').create('v1', 'Variables');
TPA3.sol('sol1').create('e1', 'Eigenvalue');
TPA3.sol('sol1').feature('e1').create('d1', 'Direct');

TPA3.result.numerical.create('int1', 'IntSurface');
TPA3.result.numerical('int1').selection.all;
TPA3.result.numerical('int1').set('probetag', 'none');
TPA3.result.create('pg1', 'PlotGroup2D');
TPA3.result('pg1').create('surf1', 'Surface');

modefreq = 'c_const/lam[um]';
modefreq = strrep(modefreq, 'lam', num2str(lam));
TPA3.study('std1').feature('mode').set('modeFreq', modefreq);
TPA3.study('std1').feature('mode').set('neigsactive', true);
TPA3.study('std1').feature('mode').set('neigs', 2);
TPA3.study('std1').feature('mode').set('shiftactive', true);
TPA3.study('std1').feature('mode').set('shift', '2.5');

TPA3.sol('sol1').attach('std1');
TPA3.sol('sol1').feature('e1').set('rtol', 1.0E-12);
TPA3.sol('sol1').feature('e1').set('transform', 'effective_mode_index');
TPA3.sol('sol1').feature('e1').set('neigs', 1);
TPA3.sol('sol1').feature('e1').set('shift', '2.5');
TPA3.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
TPA3.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (ewfd)');
TPA3.sol('sol1').runAll;

TPA3.result.numerical('int1').set('looplevelinput', {'manual'});
TPA3.result.numerical('int1').set('table', 'tbl1');
TPA3.result.numerical('int1').set('expr', {'ewfd.normE^2' 'sqrt(realdot(ewfd.Ex,ewfd.Ex)+realdot(ewfd.Ey,ewfd.Ey))^2'});
TPA3.result.numerical('int1').set('unit', {'kg^2*m^4/(s^6*A^2)' 'kg^2*m^4/(s^6*A^2)'});
TPA3.result.numerical('int1').set('descr', {'' ''});
TPA3.result.numerical('int1').setResult;
TPA3.result('pg1').label('Electric Field (ewfd)');
TPA3.result('pg1').set('looplevel', [1]);
TPA3.result('pg1').set('frametype', 'spatial');
TPA3.result('pg1').feature('surf1').set('expr', 'ewfd.Ex');
TPA3.result('pg1').feature('surf1').set('descr', 'Electric field, x component');
TPA3.result('pg1').feature('surf1').set('resolution', 'normal');

out = TPA3;
