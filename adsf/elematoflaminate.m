function [ke,me] = elematoflaminate

kmbf = guasse(@laminatembf,[2,2]);
kcc = guasse(@laminateshear,[1,1]);
ke = kmbf+kcc;
me = guasse(@laminateme,[2,2]);
