function [Like]=logitL(b,X);

Like=1 ./ (1+exp(-X*b));