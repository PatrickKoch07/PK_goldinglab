function recon_bestcv_max3proj(ip, frm, cdir, cnum, cnum_focus, crg)
%reconstruct image: 1 channel only
%Variables: ip, 
%           frame#, 
%           channel directory, 
%           channel#, 
%           LUTS (can set to [] if not sure)

%subplot 1: cell outline + raw colored image
%subplot 2: cell outline + reconstructed colored image (base on spot data)
%subplot 3: cell outline + raw colored image + spots in ellipses

%varargin: enables the function to accept any number of input arguments


%Create the cell outline
load([ip.seg.dir ip.image.base_name 'seg' num2str(frm, '%03d') '.mat'], 'LcFull');
ALL_cellPerim = zeros(size(LcFull));
for iCell = 1:1:max(LcFull(:))
    ALL_cellPerim = ALL_cellPerim + bwperim(LcFull == iCell);
end

%Create the figure
figure('numbertitle', 'off', 'name' , ['Reconstruction: frame No ' num2str(frm)], ...
    'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.95 0.75]);
ha = tight_subplot(1, 3, [.01 .01], [.1 .01], [.01 .01]);

%Remove a subplot containing raw phase image: trust my schnitzcells results
%Subplot 1: Cell outline + raw colored image
axes(ha(1)); hold on;
plotMaxPwSpt(ip, frm, cdir, cnum, 'c', 0, cnum_focus, crg);
spy(ALL_cellPerim, 'w', 0.5);
xlabel('');
axis square;
title('Raw colored image');
daspect([1 1 1]); imcontrast

%Subplot 2: Cell outline + reconstructed colored image
axes(ha(2)); hold on;
if isempty(crg)
%     MaxP = getMaxP(ip, frm, cnum);
    % get the best range
%     dMaxP = double(MaxP(MaxP<2^15)); % keep away those saturated pixels
%     bins = linspace(min(dMaxP(:)), max(dMaxP(:)), 100);
%     [ya, xa] = hist(dMaxP(:), bins);
%     maxx = setxlim(xa, ya, 0.993);
% 
%     crg2 = [mean(dMaxP(:)) maxx];
    plotReconsImg(ip, frm, cdir, cnum, LcFull, cnum_focus, crg);    %input LcFull for background determination
else
    plotReconsImg(ip, frm, cdir, cnum, LcFull, cnum_focus, crg);    %input LcFull for background determination
end
spy(ALL_cellPerim, 'w', 0.5);
xlabel('');
axis square;
title('Reconstructed fluor image from spots');
daspect([1 1 1]); imcontrast

%Subplot 3: Cell outline + raw colored image + spots in ellipses
% {
axes(ha(3)); hold on;
plotMaxPwSpt(ip, frm, cdir, cnum, 'c', 1, cnum_focus, crg);
spy(ALL_cellPerim, 'w', 0.5);
xlabel('');
axis square;
title('Ellipses showing spots');
daspect([1 1 1]); imcontrast;
%}

linkaxes(ha, 'xy');
end

function plotMaxPwSpt(ip, frm, cdir, cnum, ccol, idx, cnum_focus, varargin)
% function plot individual channel
% figure;
% imshow maxP
MaxP = getMaxP(ip, frm, cnum, cnum_focus);

% get the best range
dMaxP = double(MaxP(MaxP<2^15)); % keep away those saturated pixels
bins = linspace(min(dMaxP(:)), max(dMaxP(:)), 100);
[ya, xa] = hist(dMaxP(:), bins);
maxx = setxlim(xa, ya, 0.993);

if isempty(varargin)
    imshow(MaxP, [mean(dMaxP(:)) maxx]); hold all;
else
    imshow(MaxP, varargin{1}); hold all;
end

%show spot data only when idx == 1
if idx
    load([cdir '\data\peakdata' num2str(frm, '%03d') '.mat']) ;
    sptD = peakdata; clear peakdata;
    %cleanse the NaN, remove the complex part for now
    sptD = real(sptD);
    delete_lst = logical(sum(isnan(sptD), 2));
    sptD(delete_lst, :) = [];

%     %Oversimplified version: simple marker indicating spot
%     plot(sptD(:, 2), sptD(:, 3), '.', 'color', ccol, 'markersize', 12);
    
    %Complicated version: ellipse indicating spot area
    for ispot = 1:numel(sptD(:, 1))
        XPos = sptD(ispot, 2);
        YPos = sptD(ispot, 3);
        minor_ax = sptD(ispot, 15);
        major_ax = sptD(ispot, 16);
        %determine the angle following TY notes on 20181109:
        if sptD(ispot, 4) > sptD(ispot, 5)      % a - b > 0
            angle = 0.5 * atan(2*sptD(ispot, 6) / (sptD(ispot, 4) - sptD(ispot, 5)));
        else                                    % a - b <= 0
            angle = 0.5 * atan(2*sptD(ispot, 6) / (sptD(ispot, 4) - sptD(ispot, 5))) + pi / 2;
        end
        if pi*minor_ax*major_ax < 500
            ell = DrawEllipse(XPos, YPos, minor_ax, major_ax, angle, ccol);
        end
    end
    
    if 0
        % write spot area
        for i3=1:length(sptD(:,2))
            text((sptD(i3,2))+3*0.5,...
                (sptD(i3,3)),...
                num2str(pi*sptD(i3,16).*sptD(i3,15), '%3.1f'),...
                'Color',[155 187 89]/255,...
                'FontSize',8,'FontName','Calibri','FontWeight','Bold') ;
        end
    end
    if 0
        % write major/minor
        for i3=1:length(sptD(:,2))
            text((sptD(i3,2))+3*0.5,...
                (sptD(i3,3)),...
                num2str(sptD(i3,16)./sptD(i3,15), '%3.1f'),...
                'Color',[155 187 89]/255,...
                'FontSize',8,'FontName','Calibri','FontWeight','Bold') ;
        end
    end
end

end

function plotReconsImg(ip, frm, cdir, cnum, LcFull, cnum_focus, varargin)
%Read the fluor image for background
%Some channels has higher background signal in cells, that need to be taken into account
MaxP = getMaxP(ip, frm, cnum, cnum_focus);
in_cell_idx = logical(LcFull);          %intracellular pixels
ex_cell_idx = ~logical(LcFull);         %extracellular pixels
BG = median(double(MaxP(ex_cell_idx)));         %extracellular background
rec_img = BG * ones(ip.image.size, ip.image.size);
BG_in = median(double(MaxP(in_cell_idx)));      %intracellular background
rec_img(in_cell_idx) = BG_in;
% LUT_dev = BG_in - BG;                           %Regulate upper LUTs using this difference

load([cdir '\data\peakdata' num2str(frm, '%03d') '.mat']);
sptD = peakdata; clear peakdata;
%cleanse the NaN, remove the complex part for now
sptD = real(sptD);
delete_lst = logical(isnan(sptD(:,1)));
sptD(delete_lst, :) = [];

box_size = 50;         %Size of "fit box": distance in pixels that influenced by a spot

for ispot = 1:numel(sptD(:, 1))
    %Read the spot info
    PH = sptD(ispot, 1);
    xx = sptD(ispot, 2);
    yy = sptD(ispot, 3);
    aa = sptD(ispot, 4);
    bb = sptD(ispot, 5);
    cc = sptD(ispot, 6);
    
    %Find the influenced pixels
    x_rg = max(1, round(xx-box_size)):1:min(round(xx+box_size), ip.image.size);
    y_rg = max(1, round(yy-box_size)):1:min(round(yy+box_size), ip.image.size);
    x_mat = ones(numel(y_rg), 1) * x_rg;
    y_mat = y_rg' * ones(1, numel(x_rg));
    
    %Calculate the fluor intensity in the "fit box"
    fluo_ = PH .* exp(-(aa .* (x_mat-xx).^2 + bb .* (y_mat-yy).^2 + 2*cc .* (x_mat-xx) .* (y_mat-yy)));
    
    %Add the fluor intensity created from a spot to the reconstructed image
    rec_tmp = zeros(ip.image.size, ip.image.size);
    rec_tmp(y_rg, x_rg) = fluo_;
    rec_img = rec_img + rec_tmp;
    
end

% get the best range
dRecImg = double(rec_img(rec_img < 2^15)); % keep away those saturated pixels
bins = linspace(min(dRecImg(:)), max(dRecImg(:)), 100);
[ya, xa] = hist(dRecImg(:), bins);
maxx = setxlim(xa, ya, 0.993);

if isempty(varargin)
    imshow(rec_img, [mean(dRecImg(:)) maxx]); hold all;
else
%     LUTs = varargin{1};
%     LUTs(2) = LUTs(2) - LUT_dev;
%     imshow(rec_img, LUTs); hold all;
    imshow(rec_img, varargin{1}); hold all;
end

end

function MaxP = getMaxP(ip, n_frame, n_channel, cnum_focus)
z_range   = 1 : ip.image.zrange ;
image_path = [ip.exp.path  'images\'] ;
tmp_spl = ip.exp.frm2spl(n_frame) ;
load([ip.seg.dir ip.image.base_name 'seg' num2str(n_frame,'%03d') '.mat'] , 'LcFull');

N = max(LcFull(:));
        all_cv = zeros(numel(z_range),N);
        
        fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            %fprintf(1,['z-Slice ' num2str(i) ' out of ' num2str(numel(z_range)) '.' sprintf('\n')]);
            image_fullname  = [image_path 'sample_' num2str(tmp_spl,'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),'%02d') ];
            slice_name = [image_fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(cnum_focus, '%01d') '.tif'];
            FluorImage = imread(slice_name);
            
            for ii = 1:N
                cell = immultiply(LcFull == ii, FluorImage);
                all_cv(i, ii) = sqrt(var(double(cell(LcFull(:) == ii)))) / mean(double(cell(LcFull(:) == ii)));
            end
        end
        
        [~, max_cv_z] = max(all_cv);
        max_FluorImage = zeros(size(LcFull));
        fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            image_fullname  = [image_path 'sample_' num2str(tmp_spl,'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),'%02d') ];
            if i ~= 1
                FluorImage_past = FluorImage_curr;
                FluorImage_curr = FluorImage_fut;
                if i == numel(z_range)
                    FluorImage_fut = zeros(size(FluorImage_curr));
                else
                    slice_name = [image_fullname 'z' num2str(z_range(i+1), ...
                        '%01d') 'c' num2str(n_channel, '%01d') '.tif'];
                    FluorImage_fut = imread(slice_name);
                end
            else
                slice_name_c = [image_fullname 'z' num2str(z_range(i+1), ...
                    '%01d') 'c' num2str(n_channel, '%01d') '.tif'];
                slice_name_f = [image_fullname 'z' num2str(z_range(i+1), ...
                    '%01d') 'c' num2str(n_channel, '%01d') '.tif'];
                FluorImage_curr = imread(slice_name_c);
                FluorImage_fut = imread(slice_name_f);
                FluorImage_past = zeros(size(FluorImage_curr));
            end
            
            for ii = 1:N
                if max_cv_z(ii) == i
                    img_p = double(immultiply(LcFull == ii, FluorImage_past));
                    img_c = double(immultiply(LcFull == ii, FluorImage_curr));
                    img_f = double(immultiply(LcFull == ii, FluorImage_fut));
                    cell_proj = max(cat(3,img_p,img_c,img_f),[],3);
                    max_FluorImage = max_FluorImage + cell_proj;
                end
            end
        end
        MaxP = max_FluorImage;

end

function ell = DrawEllipse(XPos, YPos, minor_ax, major_ax, angle, color)
r = 0:0.01:2*pi;                                        %Full circle
p = [(minor_ax * cos(r))' (major_ax * sin(r))'];        %Ellipse before rotate and move
alpha = [ cos(angle) sin(angle);
         -sin(angle) cos(angle)];
p1 = p * alpha;                                         %Ellipse after rotate
p2 = p1 + [XPos YPos];                                  %Ellipse after move
ell = patch('xdata', p2(:, 1), 'ydata', p2(:, 2), 'facecolor', color, 'edgecolor', 'none');
end

