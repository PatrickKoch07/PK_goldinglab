<<<<<<< HEAD
function [CellSize,AvgCellFluor,cell_length,backgroundfluor] = FluorMeasure_indivZ(sr, Mask)

%A function measuring the fluorescence of numbered cells in an image
%
%Input:
% Fluor:: the image to be measured
% Mask:: a label matrix specifying the disjoint regions to measure
%
%Output:
% CellSize:: the area of each cell in pixels
% AvgCellFluor:: fluorescence per pixel, BG subtracted
% cell_length:: cell length

z_range = sr.image.zrange;
N = max(double(Mask(:)));

all_cv = zeros(numel(z_range),N);
fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
for i = 1:numel(z_range)
    slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
    FluorImage = imread(slice_name);

    for ii = 1:N
        cell = immultiply(Mask == ii, FluorImage);
        all_cv(i, ii) = sqrt(var(double(cell(Mask(:) == ii)))) / mean(double(cell(Mask(:) == ii)));
    end
end
[~, max_cv_z] = max(all_cv);
Fluor = zeros(size(Mask));
bgFluor = zeros(max(sr.image.zrange), 1);
fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
for i = 1:numel(z_range)
    slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
    FluorImage = imread(slice_name);          
    
    bgFluor(i) = median(FluorImage(~Mask));
    
    for ii = 1:N
        if max_cv_z(ii) == i
            Fluor = Fluor + double(immultiply(Mask == ii, FluorImage));
        end
    end
end

CellMap = Mask;
nnn=double(Mask);
N=max(nnn(:));

%Preallocate Vectors
CellSize = double(zeros(1,N));
AvgCellFluor = double(zeros(1,N));
backgroundfluor = double(zeros(1,N));
cell_length = double(zeros(1,N));
%____________________


%Loop over cells and calculate Fluorescences inside
for index=1:N
    
    Cell = immultiply( CellMap == index, Fluor);
    TotalCellFluor = sum(Cell(:));   % before substracting the background
    
    T = (CellMap == index);
    CellSize(index) = sum(T(:));
    
    backgroundfluor(index) = bgFluor(max_cv_z(N));
    
    AvgCellFluor (index) =  (TotalCellFluor - backgroundfluor(index)) / CellSize(index);
    
    % calculate cell length
    CellStats = regionprops(double(Mask==index),'Orientation','Centroid','MajorAxisLength','PixelList') ;
    alpha = -CellStats.Orientation ; % angle between x-axis and cell axis in degrees (-90 to 90) (- sign because y-axis is reversed)
    CellMAL = CellStats.MajorAxisLength ; % cell axis length from regionprops function (not so accurate)
    CellXY = CellStats.Centroid ; % cell centroid
    CellPixels = CellStats.PixelList ; % cell Pixels
    axis_line = [ CellXY(1)+[-CellMAL:1:CellMAL]*cosd(alpha) ; CellXY(2)+[-CellMAL:1:CellMAL]*sind(alpha) ] ; % line through cell centroid with same orientation as cell
    % throw away points in axis_line that are outside the cell
    for i2 = size(axis_line,2):-1:1
       D = sqrt( (axis_line(1,i2)-CellPixels(:,1)).^2 + (axis_line(2,i2)-CellPixels(:,2)).^2 ) ; % vector containing distances from all cell pixels to a point in axis_line
       if min(D)>=1
          axis_line(:,i2) = [] ;
       end
    end
    cell_length(index) = sqrt( (axis_line(1,1)-axis_line(1,end))^2 + (axis_line(2,1)-axis_line(2,end))^2 ) ; % a more accurate estimate of cell length
    
end
%__________________________________________________
=======
function [CellSize,AvgCellFluor,cell_length,backgroundfluor] = FluorMeasure_indivZ(sr, Mask)

%A function measuring the fluorescence of numbered cells in an image
%
%Input:
% Fluor:: the image to be measured
% Mask:: a label matrix specifying the disjoint regions to measure
%
%Output:
% CellSize:: the area of each cell in pixels
% AvgCellFluor:: fluorescence per pixel, BG subtracted
% cell_length:: cell length

z_range = sr.image.zrange;
N = max(double(Mask(:)));

all_cv = zeros(numel(z_range),N);
fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
for i = 1:numel(z_range)
    slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
    FluorImage = imread(slice_name);

    for ii = 1:N
        cell = immultiply(Mask == ii, FluorImage);
        all_cv(i, ii) = sqrt(var(double(cell(Mask(:) == ii)))) / mean(double(cell(Mask(:) == ii)));
    end
end
[~, max_cv_z] = max(all_cv);
Fluor = zeros(size(Mask));
bgFluor = zeros(max(sr.image.zrange), 1);
fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
for i = 1:numel(z_range)
    slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
    FluorImage = imread(slice_name);          
    
    bgFluor(i) = median(FluorImage(~Mask));
    
    for ii = 1:N
        if max_cv_z(ii) == i
            Fluor = Fluor + double(immultiply(Mask == ii, FluorImage));
        end
    end
end

CellMap = Mask;
nnn=double(Mask);
N=max(nnn(:));

%Preallocate Vectors
CellSize = double(zeros(1,N));
AvgCellFluor = double(zeros(1,N));
backgroundfluor = double(zeros(1,N));
cell_length = double(zeros(1,N));
%____________________


%Loop over cells and calculate Fluorescences inside
for index=1:N
    
    Cell = immultiply( CellMap == index, Fluor);
    TotalCellFluor = sum(Cell(:));   % before substracting the background
    
    T = (CellMap == index);
    CellSize(index) = sum(T(:));
    
    backgroundfluor(index) = bgFluor(max_cv_z(N));
    
    AvgCellFluor (index) =  (TotalCellFluor - backgroundfluor(index)) / CellSize(index);
    
    % calculate cell length
    CellStats = regionprops(double(Mask==index),'Orientation','Centroid','MajorAxisLength','PixelList') ;
    alpha = -CellStats.Orientation ; % angle between x-axis and cell axis in degrees (-90 to 90) (- sign because y-axis is reversed)
    CellMAL = CellStats.MajorAxisLength ; % cell axis length from regionprops function (not so accurate)
    CellXY = CellStats.Centroid ; % cell centroid
    CellPixels = CellStats.PixelList ; % cell Pixels
    axis_line = [ CellXY(1)+[-CellMAL:1:CellMAL]*cosd(alpha) ; CellXY(2)+[-CellMAL:1:CellMAL]*sind(alpha) ] ; % line through cell centroid with same orientation as cell
    % throw away points in axis_line that are outside the cell
    for i2 = size(axis_line,2):-1:1
       D = sqrt( (axis_line(1,i2)-CellPixels(:,1)).^2 + (axis_line(2,i2)-CellPixels(:,2)).^2 ) ; % vector containing distances from all cell pixels to a point in axis_line
       if min(D)>=1
          axis_line(:,i2) = [] ;
       end
    end
    cell_length(index) = sqrt( (axis_line(1,1)-axis_line(1,end))^2 + (axis_line(2,1)-axis_line(2,end))^2 ) ; % a more accurate estimate of cell length
    
end
%__________________________________________________
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
