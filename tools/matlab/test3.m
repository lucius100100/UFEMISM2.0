clc
clear all
close all

resid = [];
resid( end+1) =    2.0000000000000000     ;
resid( end+1) =    1.3653364942087418     ;
resid( end+1) =    1.0429772142009530     ;
resid( end+1) =   0.63088571750774791     ;
resid( end+1) =   0.30513346596550961     ;
resid( end+1) =   0.16092949054147279     ;
resid( end+1) =    9.3997348940390663E-002;
resid( end+1) =    5.8760845058691648E-002;
resid( end+1) =    2.9202198706365541E-002;
resid( end+1) =    1.7220609277156698E-002;
plot( log( resid) / log( 10))