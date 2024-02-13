%% This MCX simulation computes Jacobian for LSCI for various wavelengths
% Note that the Jacobian depends on the source type and detector type
%% intialize paths
addpath('../mcxlab-win-x86_64-v2020/mcxlab/');
addpath('../mcxlab-win-x86_64-v2020/mcx/');
addpath('../mcxlab-win-x86_64-v2020/mcx/utils');
fprintf("MCXLab code files added to path. Displaying GPU info\n"); 
mcxlab('gpuinfo');

%% specify volume size (in voxels) 
Nr = 121;
Nc = 121;
Nzr = 91; 

voxel_size_mm = 0.3;            % MCX simulation voxel size in mm
Nz = floor(12/voxel_size_mm);   % generate kernel for up to 12 mm depth 

%% scene description (scattering params)
g = 0.9;                                % scattering anisotropy factor
refr_idx =  1.37;                       % refractive index
lambda = [450, 515, 658, 850] .* 1e-6;  % in mm
mu_a = [0.04, 0.03, 0.02, 0.018];       % per mm 
mu_s = [0.95, 0.9, 0.8, 0.6] ./ (1-g);  % per mm

k_vals = 2*pi ./ lambda;                % wave number
Db = 0.5e-6;                            % diffusion coeff

%% photon replay params
t_points    = 10;
exposure_s  = 5e-4;     % 500 us
kernel      = zeros(numel(lambda), Nr, Nc, Nz);    % instantiate Jacobian
tau_vals    = linspace(0, exposure_s, t_points);

src_pos     = [(Nr-1)/2+1 (Nc-1)/2+1 0]; 
det_pos     = [(Nr-1)/2+1 (Nc-1)/2+1 1 1]; 

idx = 1;
for lmbd = lambda
    % set up configuration
    k0                  = k_vals(idx);
    cfg.seed            = 0; 
    cfg.gpuid           = 1; 
    cfg.unitinmm        = voxel_size_mm; 
    cfg.isnormalized    = 0;
    cfg.vol             = uint8(ones(Nr,Nc,Nzr)); 
    cfg.isspecular      = 0; 
    cfg.prop            = [0.0, 0.0, 1.0, 1; 
                            mu_a(idx), mu_s(idx), g, refr_idx]; 
    cfg.srcpos          = [src_pos(1) src_pos(2) src_pos(3)]; 
    cfg.srctype         = 'disk'; 
    cfg.srcparam1       = [60 50 0 360];        % radius of disk source, plus angular coverage of disk 
    cfg.srcdir          = [0 0 1];              % pointing in z-direction    
    cfg.nphoton         = 1e9; 
    cfg.isreflect       = 0; 
    cfg.issrcfrom0      = 0; 
    cfg.tstart          = 0; 
    cfg.tend            = 5e-7; 
    cfg.tstep           = 5e-7; 
    cfg.issaveexit      = 0; 
    cfg.autopilot       = 1; 
    cfg.maxdetphoton    = 3e5; 
    cfg.savedetflag     = 'dpm';
    cfg.detpos          = [det_pos(1), det_pos(2), det_pos(3), det_pos(4)]; 
    
    % run the Monte Carlo simulation
    [flux, detp, vol, seeds] = mcxlab(cfg);
    % re-weight the photons
    w1 = mcxdetweight(detp, cfg.prop); 

    % setup for photon replay 
    newcfg              = cfg; 
    newcfg.seed         = seeds.data;
    newcfg.outputtype   = 'wl'; 
    newcfg.replaydet    = 1; 
    newcfg.isnormalized = 0; 

    [g1_norm,G1]=generate_g1(mu_s(idx).*(1-g),mu_a(idx),tau_vals,'brownian',Db, lmbd, detp,voxel_size_mm);

    count = 1;
    kappa = zeros(Nr, Nc, Nz, t_points); 
    for tau=tau_vals
        newcfg.detphotons   = detp.data;
        musp                = cfg.prop(2,2) * (1 - cfg.prop(2,3)); 
        mu_tau              = 2*k0^2*musp*tau*Db; 
        muad                = cfg.prop(2,1) + mu_tau; 
        newcfg.prop(2,1)    = muad;

        % run the replayed Monte Carlo simuation for a particular tau
        [flux2, detp2, vol2, seeds2] = mcxlab(newcfg); 
        % re-weight the photons
        w2 = mcxdetweight(detp, newcfg.prop); 

        kt = (voxel_size_mm)^3 * flux2.data(:,:,1:Nz)/g1_norm; 
        kappa(:,:,:,count) = kt * G1(count) * (1-tau/exposure_s) * tau; 
        count = count+1; 
    end
    kernel(idx,:,:,:) = 1/exposure_s * sum(kappa,4) * mean(diff(tau_vals)); 
    idx = idx+1;
end
% fprintf('saving the results... \n');
% fname = sprintf('lsci_lambda_analyse_v0.mat');
% save(fname,'kernel','lambda','tau_vals','mu_a','mu_s','g','voxel_size_mm','cfg','-v7.3');
clearvars newcfg flux flux2 seeds2 detp2 seeds detp kappa w1 w2 kt cfg

fprintf('plotting the results... \n'); 
figure;
i=1; 
for idx = 1:numel(lambda)
    subplot(numel(lambda),2,i); imagesc(squeeze(log10(kernel(idx,:,:,1)))); axis image; colorbar; title("Jacobian @ z=0")
    subplot(numel(lambda),2,i+1); imagesc(squeeze(log10(kernel(idx,:,61,:)))); axis image; colorbar; title("Jacobian @ y=61"); 
    i = i+2;
end
