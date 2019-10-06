clear; close all

delta = 0:0.001:0.5;

Up = delta.^2;
Low = 2.* Up - 1;

RegPar = Up./Low;

plot(delta,RegPar)