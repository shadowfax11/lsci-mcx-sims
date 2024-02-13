addpath('../mcxlab-win-x86_64-v2020/mcxlab/');
addpath('../mcxlab-win-x86_64-v2020/mcx/');
addpath('../mcxlab-win-x86_64-v2020/mcx/utils');
fprintf("MCXLab code files added to path. Displaying GPU info\n"); 
mcxlab('gpuinfo');

Nr = 121;
Nc = 121;
Nzr = 91; 

unitinmm = 0.3;     % MCX simulation voxel size in mm
Nz = floor(12/unitinmm); % we want to generate a kernel only for <= 12mm depth

% Scene description
mu_a = 0.005;
mu_s = 7;
g = 0.9;
lamda = 660*1e-6;
k0 = 2*pi/lamda;
Db = 0.5e-6;

% Photon replay
t_points    = 10; 
EXPOSURE    = 5e-4;    % 500 us
kernel      = zeros(Nr, Nc, Nzr); 
Te          = EXPOSURE; 
TAU         = linspace(0, Te, t_points); 

for lx = 61 
    % For one source-detector configuration
    % base simulation
    cfg.seed            = 0; 
    cfg.gpuid           = 1; 
    cfg.unitinmm        = unitinmm; 
    cfg.isnormalized    = 0;
    cfg.vol             = uint8(ones(Nr,Nc,Nzr)); 
    cfg.isspecular      = 0; 
    cfg.prop            = [0.0, 0.0, 1.0, 1; 
                            mu_a, mu_s, g, 1.37]; 
    cfg.srcpos          = [lx (Nc-1)/2+1 0]; 
    cfg.srctype         = 'disk'; 
    cfg.srcparam1       = [60 50 0 360];    % radius of disk source
    cfg.srcdir          = [0 0 1];          % pointing in z-direction    
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
    cfg.detpos          = [(Nr-1)/2+1 (Nc-1)/2+1 1 1]; 
    
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

    [g1_norm,G1]=generate_g1(mu_s.*(1-0.9),mu_a,TAU,'brownian',Db, lamda, detp,unitinmm);

    count = 1;
    kappa = zeros(Nr, Nc, Nz, t_points); 
    for tau=TAU
        newcfg.detphotons   = detp.data;
        musp                = cfg.prop(2,2) * (1 - cfg.prop(2,3)); 
        mu_tau              = 2*k0^2*musp*tau*Db; 
        muad                = cfg.prop(2,1) + mu_tau; 
        newcfg.prop(2,1)    = muad;

        % run the replayed Monte Carlo simuation for a particular tau
        [flux2, detp2, vol2, seeds2] = mcxlab(newcfg); 
        % re-weight the photons
        w2 = mcxdetweight(detp, newcfg.prop); 

        kt = (unitinmm)^3 * flux2.data(:,:,1:Nz)/g1_norm; 
        kappa(:,:,:,count) = kt * G1(count) * (1-tau/Te) * tau; 
        count = count+1; 
    end
    kernel = 1/Te * sum(kappa,4) * mean(diff(TAU)); 
end
fprintf('saving the results... \n')
%     fname = sprintf('14milk_kernel_exp_%.2f_scale_%.2f_mus_%.2f_wl.mat',(EXP(T)*1e3),unitinmm,mu_s);
%     save(fname,'kernel','scene','cfg','-v7.3');
clearvars newcfg flux flux2 seeds2 detp2 seeds detp kappa w1 w2 kt cfg

figure;
subplot(1,2,1); imagesc(squeeze(log10(kernel(:,:,1)))); axis image; colorbar; title("Jacobian @ z=0")
subplot(1,2,2); imagesc(squeeze(log10(kernel(:,61,:)))); axis image; colorbar; title("Jacobian @ y=61"); 
