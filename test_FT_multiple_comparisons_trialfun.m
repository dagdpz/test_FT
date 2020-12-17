function [trl] = test_FT_multiple_comparisons_trialfun(cfg)

trl(1:cfg.trialdef.numtrl,1) = 1:200:cfg.trialdef.numtrl*200; 
trl(1:cfg.trialdef.numtrl,2) = 200:200:cfg.trialdef.numtrl*200;
trl(1:cfg.trialdef.numtrl,3) = - 100; 
trl(1:cfg.trialdef.numtrl,4) = cfg.trialdef.type*ones(1,cfg.trialdef.numtrl);


