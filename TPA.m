function out = TPA(n_Si, n_SiO2, lam)
%
% TPA.m
%
% Model exported on Nov 3 2017, 15:06 by COMSOL 5.3.0.260.

import com.comsol.model.*
import com.comsol.model.util.*

TPA = ModelUtil.create('Model');

TPA.modelPath('C:\Users\jc3e16\Documents\TPA a-Si-H waveguide');

TPA.comments(['Untitled\n\n']);

TPA.component.create('comp1', true);

TPA.component('comp1').geom.create('geom1', 2);

TPA.component('comp1').mesh.create('mesh1');

TPA.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');

TPA.study.create('std1');
TPA.study('std1').create('mode', 'ModeAnalysis');
TPA.study('std1').feature('mode').set('ngen', '5');
TPA.study('std1').feature('mode').activate('ewfd', true);

TPA.component('comp1').geom('geom1').run('');
TPA.component('comp1').geom('geom1').create('r1', 'Rectangle');
TPA.component('comp1').geom('geom1').feature('r1').set('type', 'solid');
TPA.component('comp1').geom('geom1').feature('r1').set('base', 'corner');
TPA.component('comp1').geom('geom1').feature('r1').set('pos', {'-0.7' '-0.2'});
TPA.component('comp1').geom('geom1').feature('r1').set('size', {'1.35' '0.6'});
TPA.component('comp1').geom('geom1').run('r1');
TPA.component('comp1').geom('geom1').lengthUnit('nm');
TPA.component('comp1').geom('geom1').feature('r1').set('size', [500 220]);
TPA.component('comp1').geom('geom1').feature('r1').set('pos', [-250 0]);
TPA.component('comp1').geom('geom1').run('r1');
TPA.component('comp1').geom('geom1').run('r1');
TPA.component('comp1').geom('geom1').create('r2', 'Rectangle');
TPA.component('comp1').geom('geom1').feature('r2').set('type', 'solid');
TPA.component('comp1').geom('geom1').feature('r2').set('base', 'corner');
TPA.component('comp1').geom('geom1').feature('r2').set('size', [5000 1950]);
TPA.component('comp1').geom('geom1').feature('r2').set('pos', [-2500 -1950]);
TPA.component('comp1').geom('geom1').runPre('fin');
TPA.component('comp1').geom('geom1').run('r2');
TPA.component('comp1').geom('geom1').create('r3', 'Rectangle');
TPA.component('comp1').geom('geom1').feature('r3').set('type', 'solid');
TPA.component('comp1').geom('geom1').feature('r3').set('base', 'corner');
TPA.component('comp1').geom('geom1').feature('r3').set('pos', {'-6000' '-4000'});
TPA.component('comp1').geom('geom1').feature('r3').set('size', {'12000' '8000'});
TPA.component('comp1').geom('geom1').run('r3');
TPA.component('comp1').geom('geom1').runPre('fin');
TPA.component('comp1').geom('geom1').feature('r1').label('Waveguide');
TPA.component('comp1').geom('geom1').feature('r2').label('Substrate');
TPA.component('comp1').geom('geom1').feature('r3').label('Air');
TPA.component('comp1').geom('geom1').run;

TPA.component('comp1').material.create('mat1', 'Common');
TPA.component('comp1').material('mat1').label('Silicon');
TPA.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
TPA.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'3.5'});
TPA.component('comp1').material('mat1').selection.set([3]);
TPA.component('comp1').material.create('mat2', 'Common');
TPA.component('comp1').material('mat2').label('Silica');
TPA.component('comp1').material('mat2').selection.set([2]);
TPA.component('comp1').material.create('mat3', 'Common');
TPA.component('comp1').material('mat3').label('Air');
TPA.component('comp1').material('mat3').selection.set([1]);
TPA.component('comp1').material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');
TPA.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', {'1'});
TPA.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
TPA.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'1.444'});

modefreq = 'c_const/lam[um]';
modefreq = strrep(modefreq, 'lam', num2str(lam));
TPA.study('std1').feature('mode').set('modeFreq', modefreq);

TPA.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', n_Si);
TPA.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', n_SiO2);

TPA.component('comp1').mesh('mesh1').autoMeshSize(1);
TPA.component('comp1').mesh('mesh1').run;

TPA.study('std1').feature('mode').set('shiftactive', true);
TPA.study('std1').feature('mode').set('shift', '3.5');
TPA.study('std1').feature('mode').set('neigsactive', true);
TPA.study('std1').feature('mode').set('neigs', 1);

TPA.sol.create('sol1');
TPA.sol('sol1').study('std1');

TPA.study('std1').feature('mode').set('notlistsolnum', 1);
TPA.study('std1').feature('mode').set('notsolnum', '1');
TPA.study('std1').feature('mode').set('listsolnum', 1);
TPA.study('std1').feature('mode').set('solnum', '1');

TPA.sol('sol1').create('st1', 'StudyStep');
TPA.sol('sol1').feature('st1').set('study', 'std1');
TPA.sol('sol1').feature('st1').set('studystep', 'mode');
TPA.sol('sol1').create('v1', 'Variables');
TPA.sol('sol1').feature('v1').set('control', 'mode');
TPA.sol('sol1').create('e1', 'Eigenvalue');
TPA.sol('sol1').feature('e1').set('neigs', 6);
TPA.sol('sol1').feature('e1').set('shift', '1');
TPA.sol('sol1').feature('e1').set('control', 'mode');
TPA.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
TPA.sol('sol1').feature('e1').create('d1', 'Direct');
TPA.sol('sol1').feature('e1').feature('d1').set('linsolver', 'mumps');
TPA.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (ewfd)');
TPA.sol('sol1').attach('std1');

TPA.result.create('pg1', 'PlotGroup2D');
TPA.result('pg1').label('Electric Field (ewfd)');
TPA.result('pg1').set('frametype', 'spatial');
TPA.result('pg1').set('data', 'dset1');
TPA.result('pg1').feature.create('surf1', 'Surface');
TPA.result('pg1').feature('surf1').set('data', 'parent');

TPA.sol('sol1').runAll;

TPA.result('pg1').run;

TPA.sol('sol1').feature('e1').set('rtol', 1.0E-12);
TPA.sol('sol1').runAll;

TPA.result('pg1').run;

out = TPA;
