function reconstructimage_v7_w_phase(ip, n_frame, n_channel, gated_folder)
%reconstructimage_v4 is function that reconstruct the image with the
%segmentation and spot recognition results. 
% main difference from previous version
% not plotting the overlapping region between gene and nasent as yellow,
% just first classify the RNA spot into nascent / cytoplasmic then plot
% nascent as yellow and cytoplasmic as red
% 
% update 022315
% now add another round of finding neighbour nascent mRNA 
% and in this round, there is a intensity-dependent distance determination
% 
% SYNOPSIS: reconstructimage_v4(ip, n_frame)
% 
% INPUT ip: the output of InitializeExp       
%       n_frame: the number of the frame to check
% OUTPUT figure
% 
% REMARKS
% 
% EXAMPLES:
% 
% 
% created with MATLAB ver.: 7.13.0.564 (R2011b) on Windows 
% 
% created by: M.W.
% DATE: 01/29/15
%
% 
%
%____________________________________________%
%% read peakdata.mat
%____________________________________________%

% load the gated spot file 
load([gated_folder '\data\' 'peakdata' num2str(n_frame, '%03d') '.mat']) ;
c3_peakdata = peakdata ;
clear peakdata    

% load the cell outline
load([ip.seg.dir ip.image.base_name 'seg' num2str(n_frame,'%03d') '.mat'] , 'LcFull');
ALL_cellPerim = zeros(size(LcFull)) ;
for iCell = 1:1:max(LcFull(:))
    ALL_cellPerim = ALL_cellPerim + bwperim(LcFull==iCell) ;
end


% plot again and mark the RNA number
% the maximum intensity projection of fluorescence channel
z_range   = 1 : ip.image.zrange ;
MaxP_c2      = zeros(ip.image.size, ip.image.size) ;

image_path = [ip.exp.path  'images\'] ;
tmp_spl = ip.exp.frm2spl(n_frame) ;


% load fluorescent images
% channel_num = n_channel;
% for i1 = 1:numel(z_range)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     image_fullname  = [image_path 'sample_' num2str(tmp_spl,'%03d') ...
%         'xy' num2str(ip.exp.frm2img(n_frame),'%02d') ];
%     slice_name = [image_fullname 'z' num2str(z_range(i1), '%01d') 'c' num2str(channel_num, '%01d') '.tif'];
%     if ~isfile(slice_name)
%         image_fullname  = [image_path 'sample_' num2str(tmp_spl,'%03d') ...
%             'xy' num2str(ip.exp.frm2img(n_frame),'%01d') ];
%         slice_name = [image_fullname 'z' num2str(z_range(i1), '%01d') 'c' num2str(channel_num, '%01d') '.tif'];
%     end
%     FluorImage = imread(slice_name);
%     MaxP_c2 = max( cat(3, MaxP_c2, FluorImage), [], 3);
% end
        N = max(LcFull(:));
        all_cv = zeros(numel(z_range),N);
        
        fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            %fprintf(1,['z-Slice ' num2str(i) ' out of ' num2str(numel(z_range)) '.' sprintf('\n')]);
            image_fullname  = [image_path 'sample_' num2str(tmp_spl,'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),'%02d') ];
            slice_name = [image_fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(1, '%01d') '.tif'];
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
            slice_name = [image_fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(n_channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);          
            
            for ii = 1:N
                if max_cv_z(ii) == i
                max_FluorImage = max_FluorImage + double(immultiply(LcFull == ii, FluorImage));
                end
            end
        end

figure('name' , ['frame No' num2str(n_frame)] , ...
       'unit' , 'normalized' , ...
       'position' , [0  0   1  1]) ;
%hold all; 

% fluoresence image
%imshow((MaxP_c2), [min(MaxP_c2(:)) max(MaxP_c2(:))]); hold on;
% imshow((MaxP_c2), [250 1500]); hold on ; 
imshow(max_FluorImage); hold on;
% cell outline
spy(ALL_cellPerim,'w',0.5); hold on ; 

% spot indexing
plot(c3_peakdata(:,2), c3_peakdata(:,3), 'ro') ;  hold on ;

xlabel('') % get rid of the "nz=...."
axis square
imcontrast
daspect([1 1 1])