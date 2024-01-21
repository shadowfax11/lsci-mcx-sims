addpath('../mcxlab-win-x86_64-v2020/mcxlab/');
addpath('../mcxlab-win-x86_64-v2020/mcx/');
addpath('../mcxlab-win-x86_64-v2020/mcx/utils');
fprintf("MCXLab code files added to path. Displaying GPU info\n"); 
mcxlab('gpuinfo');

num_Tds = 121; 
Nr = num_Tds; 
Nc = Nr; 

unitinmm = 0.3;     % MCX simulation voxel size
Nz = floor(8/unitinmm); 
Nzr = 91;

%------------Scene description------------%
src_slit = 1;
diffuse = 0;
N = Nr;
mu_a = 0.005;
mu_s = 7;
g = 0.9;
lamda = 660*1e-6;
k0 = 2*pi/lamda;
Db = 0.5e-6;

%------------Absorption reference signal------------%
min_lx = -(num_Tds-Nr)/2+1;
max_lx = Nr + (num_Tds-Nr)/2;
Lx = linspace(1, num_Tds, Nr); 
% starting illumination line is at row 50, can start from 0, but will be assymmetric

%-------------Photon replay-------------%
K = 1; 
EXPOSURE = 5e-4;    % 500 us 
t_points = 10; 

% kernel is defined for a 1x121x121x121x26 volume (???)
k_kernel = zeros(K,Nr,Nr,Nc,Nz); 
Te = EXPOSURE; 
TAU = linspace(0,Te,t_points); 
for kk=1:K
    src_count = 1; 
    for lx=40
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
        cfg.srcpos          = [lx 1 0]; 
        cfg.srctype         = 'slit'; 
        cfg.srcparam1       = [0 Nc-1 0 0];     % width of slit source
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
        k_kernel(kk,lx,:,:,:) = 1/Te * sum(kappa,4) * mean(diff(TAU)); 
        src_count = src_count + 1;
    end
    kernel = squeeze(mean(k_kernel(1:kk,:,:,:,:),1)); 
    kernel = permute(kernel, [1 4 2 3]); 
    fprintf('saving the results... \n')
%     fname = sprintf('14milk_kernel_exp_%.2f_scale_%.2f_mus_%.2f_wl.mat',(EXP(T)*1e3),unitinmm,mu_s);
%     save(fname,'kernel','scene','cfg','-v7.3');
    clearvars newcfg flux flux2 seeds2 detp2 seeds detp kappa w1 w2 kt cfg
end