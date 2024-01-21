%% Testing out MCXLAB MATLAB simulations

addpath('../mcxlab-win-x86_64-v2020/mcxlab/');
fprintf("MCXLab code files added to path. Displaying GPU info\n"); 
mcxlab('gpuinfo'); 

% define simulation
cfg.nphoton = 1e8;
cfg.vol     = ones(60,60,60,'uint8'); 
cfg.tstart  = 0; 
cfg.tend    = 5e-9; 
cfg.tstep   = 5e-9; 
cfg.srcpos  = [31,31,0]; 
cfg.srcdir  = [0,0,1]; 
cfg.prop    = [[0,0,1,1];[0.005,0.1,0.01,1.37]]; 

[flux, detp, vol, seeds]=mcxlab(cfg);
