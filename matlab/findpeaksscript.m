% FINDPEAKSSCRIPT Code associated with the LMRG Study 3
%
% Analyzes TIFF files paired with INI files as defined in PSFj.
% Identifies beads, fits a 2D Gaussian and outputs all parameters.
%
% Usage : Change the 'indir' and 'outdir' parameters and run. 
%
% Requirements : 
%    - tested with MATLAB2017b
%    - 'fullGaussian.m' and 'Gaussian1D.m'
%    - 'bfmatlab'
%
% Author : Thomas Pengo (tp81), tpengo@umn.edu, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETERS

% Directory with a list of TIFF stacks to analyze
indir = '/data/tpengo/ABRF2018/Anonymized/';

% Output directory for CSV files
outdir = '/data/tpengo/ABRF2018/Reproducible_20220728_2231CT_wN_wSNR';

% Directory to Bioformats for MATLAB
bfmatlabdir = '/home/umii/tpengo/SW/bfmatlab';

% Local window size for PSF analysis
MdataSize=[18 18 16];

%% MAIN CODE
if ~exist(outdir,'dir')
    mkdir(outdir)
end

files=dir(fullfile(indir,'*.tiff'));

import java.io.FileReader
import java.util.Properties
pxxs = zeros(length(files),1);
pzzs = zeros(length(files),1);
ris = zeros(length(files),1);
beadsizes = zeros(length(files),1);
wavelengths = zeros(length(files),1);
nas = zeros(length(files),1);
for ifile = 1:length(files)
    fname = files(ifile).name;
       
    infile = fullfile(indir,fname);

    vp = Properties();
    
    try
        vp.load(FileReader(strrep(infile,'.tiff','.ini')));

        pxxs(ifile) = str2double(vp.get('width'));
        pzzs(ifile) = str2double(vp.get('depth'));
        ris(ifile) = str2double(vp.get('ri'));
        beadsizes(ifile) = str2double(vp.get('beadsize'));
        wavelengths(ifile) = str2double(vp.get('wavelength'));
        nas(ifile) = str2double(vp.get('na'));
    catch err
        fprintf('Error opening INI file, so skipping\n');
    end
end

parfor ifile = 1:length(files)    
    addpath(bfmatlabdir)
    
    fname = files(ifile).name;
    
    pxx = pxxs(ifile);
    pzz = pzzs(ifile);
    
    if pxx>0
    

        if fname(1)~='.' && ~files(ifile).isdir
            outfile = fullfile(outdir,strrep(fname,'.tiff','.csv'));
            
            if exist(outfile,'file')
                fprintf('Skipping %s as already analyzed\n',outfile);
            else
                fprintf('Analyzing... %s\n',fname);

                infile = fullfile(indir,fname);

                v = bfopen(infile);
                v = cat(3, v{1}{:, 1});
                
                if ndims(v)<3
                    fprintf('Skipping... 2D image');
                    continue
                end

                fprintf('Peak finding..');
                
                % Create 1D bionmial kernel
                k = [1 2 1]/4;
                k = conv(k,k);
                k = conv(k,k);

                % 3 1-D convolutions along each dimension
                v1=convn(v,k,'same');
                v1=convn(v1,k','same');
                v1=convn(v1,reshape(k,1,1,length(k)),'same');

                % Find regional maxima in 3-times downsampled image
                bw = imregionalmax(v1(1:3:end,1:3:end,:));
                p = regionprops(bw,v(1:3:end,1:3:end,:),{'Area','BoundingBox','MeanIntensity','Centroid'});

                % Find mode of the image as the maximum of the histogram
                [h, hb] = hist(v1(:),100);
                [mv, imv] = max(h);
                bg = hb(imv);
                
                % Define threshold as 3 times the mode
                t = bg*3;

                % Find regional maxima whose intensity is above the
                % threshold
                mi = cat(1,p.MeanIntensity);
                q = p(mi>t);
                fprintf('done. Found %d regional maxima, filtered down to %d.\n', length(p), length(q));

                mi = cat(1,q.MeanIntensity);
                ci = cat(1,q.Centroid);

                % Start peak analysis
                fprintf('Analyzing each peak..\n');
                warning('off');
                
                fits = table;
                opts = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt');
                for i=1:length(q) 
                    
                    % Find current centroid's coordinates in full-resolution image
                    c = round(ci(i,:).*[3 3 1])
                    
                    % Find bounding box for local window around centroid
                    cmin = max(1,c-MdataSize);
                    cmax = min([size(v,2) size(v,1) size(v,3)],c+MdataSize);

                    % Get image intensities around centroid
                    ZZ = double(v(cmin(2):cmax(2),cmin(1):cmax(1),cmin(3):cmax(3)));
                    P = zeros(size(ZZ,3),1);
                    for iz=1:size(ZZ,3)
                        % Extract iz's slice
                        sl = ZZ(:,:,iz);
                        
                        % Calculate STD
                        P(iz) = std(sl(:));
                    end
                    
                    % Estimate focal plane as the slice with max STD
                    [~, imx] = max(P);
                    Z = ZZ(:,:,imx);
                    
                    % Get peak position in XY using the estimated focal
                    % plane
                    [~, I] = max(Z(:));
                    [pi, pj] = ind2sub(size(Z),I);
                    
                    % Get the Z profile along the estimated peak position
                    ZP = squeeze(double(ZZ(pi,pj,:)));
                    
                    % 1D Gaussian fit along Z
                    % [Ib, I0, z0, sz]
                    x0 = [bg, q(i).MeanIntensity-bg, 0, 5]
                    xdata = (cmin(3):cmax(3))'-imx-cmin(3);
                    [x,resnorm,residual,exitflag] = lsqcurvefit(@Gaussian1D,x0,xdata,ZP,[],[],opts);
                    sigmaz = abs(x(end));
                    
                    % FWHMz from sigma_z
                    FWHMz = 2*sqrt(2*log(2))*sigmaz;
                    
                    % Readjust plane using fitted position
                    imx1 = floor(imx+x(3));
                    if imx1>size(ZZ,3) || imx1<1
                        continue
                    end
                    Z = ZZ(:,:,imx1);
                    
                    % Plot both
                    % figure; subplot(1,2,1); imagesc(ZZ(:,:,imx)); subplot(1,2,2); imagesc(ZZ(:,:,imx1))
                    
                    % Code to plot the fit
                    % xdataf = linspace(xdata(1),xdata(end),100);ZPf = Gaussian1D(x,xdataf);figure; plot(xdata,ZP,'*k');hold;plot(xdataf,ZPf,'r'); plot(linspace(x(3)-FWHMz/2,x(3)+FWHMz/2,100),repmat(x(1)+x(2)/2,100),'-b')

                    % Perform fit on estimated focal plane
                    x = cmin(1):cmax(1); x = x-c(1);
                    y = cmin(2):cmax(2); y = y-c(2);
                    [X,Y] = meshgrid(x,y);
                    xdata = zeros(size(X,1),size(Y,2),2);
                    xdata(:,:,1) = X;
                    xdata(:,:,2) = Y;

                    x0 = [bg, q(i).MeanIntensity-bg, 0, 5, 0, 5, 0];
                    [x,resnorm,residual,exitflag] = lsqcurvefit(@fullGaussian,x0,xdata,Z,[],[],opts);

                    E = fullGaussian(x,xdata);
                    O = Z;
                    X2 = sum((O(:)-E(:)).^2./E(:));

                    Zf = fullGaussian(x,xdata);

                    P = corrcoef([Z(:) Zf(:)]);
                    
                    top_half = (Zf-x(1))>=(x(2)/2);
                    signal_top_half = mean(Z(top_half)-x(1));
                    noise_top_half = std(Z(top_half)-Zf(top_half));
                    snr_top_half = signal_top_half / noise_top_half;
                    n_top_half = sum(top_half(:));

                    fits{i,'Bkg'}= x(1);
                    fits{i,'Amp'}= x(2);
                    fits{i,'dX0'} = x(3)*pxx;
                    fits{i,'X0'} = (c(1)+x(3))*pxx;
                    fits{i,'sX'} = abs(x(4)*pxx);
                    fits{i,'dY0'} = x(5)*pxx;
                    fits{i,'Y0'} = (c(2)+x(5))*pxx;
                    fits{i,'sY'} = abs(x(6)*pxx);
                    fits{i,'Z'} = (imx1-1+cmin(3))*pzz;
                    fits{i,'sZ'} = sigmaz*pzz;
                    fits{i,'theta'} = mod(x(7),2*pi)*180/pi;
                    fits{i,'resnorm'} = resnorm;
                    fits{i,'R2'} = P(2);
                    fits{i,'N'} = numel(Z);
                    fits{i,'X2r'} = X2/numel(x0);
                    fits{i,'X2pvalue'} = chi2cdf(X2,numel(Z)-numel(x0),'upper');
                    fits{i,'N_top_half'} = n_top_half;
                    fits{i,'signal_top_half'} = signal_top_half;
                    fits{i,'noise_top_half'} = noise_top_half;
                    fits{i,'SNR_top_half'} = snr_top_half;
                    fits{i,'SNR_top_half_dB'} = 20*log10(snr_top_half);
                    

                    if mod(i,100)==99
                        fprintf('%s %d/%d\n',fname,i,length(q));
                    end
                end
                fprintf('\n');
                warning('on');
                
                if length(q)>1
                    fits.FWHMx = 2*sqrt(2*log(2))*fits.sX;
                    fits.FWHMy = 2*sqrt(2*log(2))*fits.sY;
                    fits.FWHMz = 2*sqrt(2*log(2))*fits.sZ;
                    fits.FWHMmin = min(fits.FWHMx,fits.FWHMy);
                    fits.FWHMmax = max(fits.FWHMx,fits.FWHMy);
                    fits.Asymmetry = fits.FWHMmax./fits.FWHMmin;

                    fits.valid = abs(fits.dX0./fits.FWHMx)<1 & abs(fits.dY0./fits.FWHMy)<1 & fits.R2>.9;

                    writetable(fits,outfile);
                end
            end
        end
    else
        fprintf('Skipping %s as invalid pixel size...',fname,pxx);
    end
end

%end
