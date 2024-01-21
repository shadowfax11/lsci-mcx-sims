function [g1_norm,g1]=generate_g1(musp,mua,tau, disp_model, DV, lambda, detpt,unitinmm)
%
%function [tau,g1]=generate_g1(fhist,tau, disp_model, DV, lambda, format)
%
%   Compute simulated electric-field auto-correlation function using
%   simulated photon pathlengths and scattering momentum transfer
%
%   author: Stefan Carp (carp <at> nmr.mgh.harvard.edu)
%
%   input:
%       fhist:      the file name of the output .mch file 
%       tau:        correlation times at which to compute g1 
%                   (default: 1e-7 to 1e-1 seconds, log equidistant)
%       disp_model: displacement model ('brownian', 'random_flow', <custom>)
%                   (default: brownian, see further explanation below)
%       disp_var:   value of displacement variable using mm as unit of
%                   length and s as unit of time
%                   (default: 1e-7 mm^2/s, see further explanation below)
%       lambda:     wavelenght of light used in nm
%                   (default: 785)
%       format:     the format used to save the .mch file 
%                   (default: 'float')
%
%   output:
%
%       tau:        correlation times at which g1 was computed provided for
%                   convenience (copied from input if set, otherwise 
%                   outputs default)
%       g1:         field auto-correlation curves, one for each detector
%
%   The displacement model indicates the formula used to compute the root
%   mean square displacement of scattering particles during a given delay
%   
%   brownian:       RMS= 6 * DV * tau; 
%                   DV(displacement variable)=Db (brownian diffusion coeff)
%   random_flow:    RMS= DV^2 * tau^2; 
%                   DV = V (first moment of velocity distribution)
%   <custom>:       any string other than 'brownian' or 'random_flow' will
%                   be evaluate as is using Matlab evalf, make sure it uses
%                   'DV' as the flow related independent variable, tau is
%                   indexed as tau(J). Any additional parameters can be 
%                   sent via "varargin"
%
%   This file is part of Mesh-Based Monte Carlo
%   License: GPLv3, see http://mcx.sf.net for details
%

detnum = max(detpt.detid);

if strcmp(disp_model,'brownian'),
    disp_str='rmsdisp=6*DV.*tau(J);';
elseif strcmp(disp_model,'random_flow'),
    disp_str='rmsdisp=DV.^2.*tau(J).^2;';
else 
    disp_str=['rmsdisp=' disp_model ';'];
end

n=1;
k0=2*pi/(lambda);

g1=zeros(detnum,length(tau));

for I=1:detnum,
    idx= find(detpt.detid(:)==I);     
    fprintf('Processing detector %.0f: %.0f photons\n',I,length(idx));
    pp = detpt.ppath(idx,:)*unitinmm;
    for J=1:length(tau),
        rmsdisp = 6*DV*tau(J);
        g1(I,J)=sum(exp(sum(-(k0.^2*rmsdisp/3).*musp.*pp,2)).*exp(-sum(mua.*pp,2)));
    end
    g1_norm=sum(exp(-sum(mua.*pp,2)));
    g1(I,:)=g1(I,:)./g1_norm;
end
    


