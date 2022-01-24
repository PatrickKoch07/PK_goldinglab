function [CellSize,AvgCellFluor,cell_length, BackgroundFluor] = FluorMeasure1(Fluor,Mask)

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

FBackground = immultiply(imerode(Mask == 0,strel('disk',10)), Fluor);

% calculate Background based on median
b = double(FBackground);
BckFluor_median = median(nonzeros(b));

% take median background
BackgroundFluor = BckFluor_median;

CellMap = Mask;
nnn=double(Mask);
N=max(nnn(:));

%Preallocate Vectors
CellSize = double(zeros(1,N));
AvgCellFluor = double(zeros(1,N));
cell_length = double(zeros(1,N));
%____________________


%Loop over cells and calculate Fluorescences inside
for index=1:N
    
    Cell = immultiply( CellMap == index, Fluor);
    TotalCellFluor = sum(Cell(:));   % before substracting the background
    
    T = (CellMap == index);
    CellSize(index) = sum(T(:));
    
    AvgCellFluor (index) =  TotalCellFluor / CellSize(index) - BackgroundFluor;
    
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
