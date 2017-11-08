function out = TPA2(n_Si, n_SiO2, lam)
% TPA2.m
%
% Model exported on Nov 7 2017, 11:15 by COMSOL 5.3.0.260.

import com.comsol.model.*
import com.comsol.model.util.*

TPA2 = ModelUtil.create('Model');

TPA2.modelPath('C:\Users\jc3e16\Documents\TPA a-Si-H waveguide');

TPA2.label('a-Si-H waveguide test 2.mph');

TPA2.comments(['Untitled\n\n']);

TPA2.component.create('comp1', true);

TPA2.component('comp1').geom.create('geom1', 2);

TPA2.result.table.create('evl2', 'Table');
TPA2.result.table.create('tbl1', 'Table');

TPA2.component('comp1').mesh.create('mesh1');

TPA2.component('comp1').geom('geom1').lengthUnit('nm');
TPA2.component('comp1').geom('geom1').scaleUnitValue(true);
TPA2.component('comp1').geom('geom1').create('r1', 'Rectangle');
TPA2.component('comp1').geom('geom1').feature('r1').label('Waveguide');
TPA2.component('comp1').geom('geom1').feature('r1').set('pos', [-250 0]);
TPA2.component('comp1').geom('geom1').feature('r1').set('size', {'500.0' '200'});
TPA2.component('comp1').geom('geom1').create('r2', 'Rectangle');
TPA2.component('comp1').geom('geom1').feature('r2').label('Substrate');
TPA2.component('comp1').geom('geom1').feature('r2').set('pos', [-5000 -1950]);
TPA2.component('comp1').geom('geom1').feature('r2').set('size', {'10000' '1950.0'});
TPA2.component('comp1').geom('geom1').create('r3', 'Rectangle');
TPA2.component('comp1').geom('geom1').feature('r3').label('Air');
TPA2.component('comp1').geom('geom1').feature('r3').set('pos', [-5000 0]);
TPA2.component('comp1').geom('geom1').feature('r3').set('size', {'10000.0' '9000.0'});
TPA2.component('comp1').geom('geom1').create('r4', 'Rectangle');
TPA2.component('comp1').geom('geom1').feature('r4').label('PML (air)');
TPA2.component('comp1').geom('geom1').feature('r4').set('pos', [-5500 0]);
TPA2.component('comp1').geom('geom1').feature('r4').set('size', {'11000.0' '9500'});
TPA2.component('comp1').geom('geom1').create('r5', 'Rectangle');
TPA2.component('comp1').geom('geom1').feature('r5').active(false);
TPA2.component('comp1').geom('geom1').feature('r5').label('Silicon base');
TPA2.component('comp1').geom('geom1').feature('r5').set('pos', [-5000 -4000]);
TPA2.component('comp1').geom('geom1').feature('r5').set('size', [10000 2050]);
TPA2.component('comp1').geom('geom1').create('r6', 'Rectangle');
TPA2.component('comp1').geom('geom1').feature('r6').active(false);
TPA2.component('comp1').geom('geom1').feature('r6').label('PML (silicon base)');
TPA2.component('comp1').geom('geom1').feature('r6').set('pos', [-5500 -4500]);
TPA2.component('comp1').geom('geom1').feature('r6').set('size', [11000 2550]);
TPA2.component('comp1').geom('geom1').create('r7', 'Rectangle');
TPA2.component('comp1').geom('geom1').feature('r7').label('PML (substrate)');
TPA2.component('comp1').geom('geom1').feature('r7').set('pos', [-5500 -2450]);
TPA2.component('comp1').geom('geom1').feature('r7').set('size', [11000 2450]);
TPA2.component('comp1').geom('geom1').run;

TPA2.component('comp1').material.create('mat1', 'Common');
TPA2.component('comp1').material.create('mat2', 'Common');
TPA2.component('comp1').material.create('mat3', 'Common');
TPA2.component('comp1').material('mat1').selection.set([5]);
TPA2.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
TPA2.component('comp1').material('mat2').selection.set([1 3]);
TPA2.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
TPA2.component('comp1').material('mat3').selection.set([2 4]);
TPA2.component('comp1').material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');

TPA2.component('comp1').coordSystem.create('pml1', 'PML');
TPA2.component('comp1').coordSystem('pml1').selection.set([1 2]);

TPA2.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
TPA2.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 1);
TPA2.component('comp1').physics('ewfd').feature('sctr1').selection.set([8 10 17]);

TPA2.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
TPA2.component('comp1').mesh('mesh1').create('ref1', 'Refine');

TPA2.result.table('evl2').label('Evaluation 2D');
TPA2.result.table('evl2').comments('Interactive 2D values');
TPA2.result.table('tbl1').comments('Surface Integration 1 (ewfd.normE^2)');

TPA2.component('comp1').view('view1').axis.set('xmin', -9760.162109375);
TPA2.component('comp1').view('view1').axis.set('xmax', 9760.162109375);
TPA2.component('comp1').view('view1').axis.set('ymin', -3047.5);
TPA2.component('comp1').view('view1').axis.set('ymax', 10097.5);
TPA2.component('comp1').view('view1').axis.set('abstractviewlratio', -0.3872874677181244);
TPA2.component('comp1').view('view1').axis.set('abstractviewrratio', 0.3872874677181244);
TPA2.component('comp1').view('view1').axis.set('abstractviewbratio', -0.05000000074505806);
TPA2.component('comp1').view('view1').axis.set('abstractviewtratio', 0.05000000074505806);
TPA2.component('comp1').view('view1').axis.set('abstractviewxscale', 21.90833282470703);
TPA2.component('comp1').view('view1').axis.set('abstractviewyscale', 21.90833282470703);

TPA2.component('comp1').material('mat1').label('Silicon');
TPA2.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
TPA2.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
TPA2.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', n_Si);
TPA2.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
TPA2.component('comp1').material('mat2').label('Silica');
TPA2.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', '');
TPA2.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('ki', '');
TPA2.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', n_SiO2);
TPA2.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
TPA2.component('comp1').material('mat3').label('Air');
TPA2.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', '');
TPA2.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', '');
TPA2.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
TPA2.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});

TPA2.component('comp1').coordSystem('pml1').set('ScalingType', 'userDefined');
TPA2.component('comp1').coordSystem('pml1').set('directions', '2');
TPA2.component('comp1').coordSystem('pml1').set('dmax', {'1e-5' '1e-5' '1'});

TPA2.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
TPA2.component('comp1').mesh('mesh1').feature('ref1').set('boxcoord', true);
TPA2.component('comp1').mesh('mesh1').feature('ref1').set('xmax', 963);
TPA2.component('comp1').mesh('mesh1').feature('ref1').set('xmin', -885);
TPA2.component('comp1').mesh('mesh1').feature('ref1').set('ymax', 601);
TPA2.component('comp1').mesh('mesh1').feature('ref1').set('ymin', -580);
TPA2.component('comp1').mesh('mesh1').run;

TPA2.study.create('std1');
TPA2.study('std1').create('mode', 'ModeAnalysis');

TPA2.sol.create('sol1');
TPA2.sol('sol1').study('std1');
TPA2.sol('sol1').attach('std1');
TPA2.sol('sol1').create('st1', 'StudyStep');
TPA2.sol('sol1').create('v1', 'Variables');
TPA2.sol('sol1').create('e1', 'Eigenvalue');
TPA2.sol('sol1').feature('e1').create('d1', 'Direct');

TPA2.result.numerical.create('int1', 'IntSurface');
TPA2.result.numerical('int1').selection.set([3 5]);
TPA2.result.numerical('int1').set('probetag', 'none');
TPA2.result.create('pg1', 'PlotGroup2D');
TPA2.result('pg1').create('surf1', 'Surface');

modefreq = 'c_const/lam[um]';
modefreq = strrep(modefreq, 'lam', num2str(lam));
TPA2.study('std1').feature('mode').set('modeFreq', modefreq);
TPA2.study('std1').feature('mode').set('neigsactive', true);
TPA2.study('std1').feature('mode').set('neigs', 1);
TPA2.study('std1').feature('mode').set('shiftactive', true);
TPA2.study('std1').feature('mode').set('shift', '3.6');

TPA2.sol('sol1').attach('std1');
TPA2.sol('sol1').feature('e1').set('rtol', 1.0E-12);
TPA2.sol('sol1').feature('e1').set('transform', 'effective_mode_index');
TPA2.sol('sol1').feature('e1').set('neigs', 1);
TPA2.sol('sol1').feature('e1').set('shift', '3.6');
TPA2.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
TPA2.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (ewfd)');
TPA2.sol('sol1').runAll;

TPA2.result.numerical('int1').set('looplevelinput', {'manual'});
TPA2.result.numerical('int1').set('table', 'tbl1');
TPA2.result.numerical('int1').set('expr', {'ewfd.normE^4'});
TPA2.result.numerical('int1').set('unit', {'kg^4*m^6/(s^12*A^4)'});
TPA2.result.numerical('int1').set('descr', {''});
TPA2.result.numerical('int1').setResult;
TPA2.result('pg1').label('Electric Field (ewfd)');
TPA2.result('pg1').set('looplevel', [1]);
TPA2.result('pg1').set('frametype', 'spatial');
TPA2.result('pg1').feature('surf1').set('expr', 'ewfd.Ex');
TPA2.result('pg1').feature('surf1').set('descr', 'Electric field, x component');
TPA2.result('pg1').feature('surf1').set('resolution', 'normal');

out = TPA2;
