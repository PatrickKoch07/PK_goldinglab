function outim = imshowlabel(L,bwscreen)
% IMSHOWLABEL is used to display an integer image 

if islogical(L),
    L = double(L);
end;

% % L2 has every non-background blob in range [2,256] and 
% % sets background to one, corresp. to first entry in mymap
% L2 = mod(L,255)+2;
% L2(L==0) = 1;

%pk (11/26/2021)
L2 = mod(L,254)+3;
L2(L==0) = 1;
L2(L<0) = 2;
% L(L < 0) = 1;

% M is the maximum color table entry, at most 256 colors
M = min(max2(L)+2,256);
% % % % create a color map
% % % mymap = hsv(M);
% % % % explicitly set the colormap's first entry to black for background
% % % mymap(1,:)=[0 0 0];

% -------------------
% begin Tommy, Mar 14, 2009
% create a color map
A = [79 129 189  ; ...
     192 80 77   ; ...
     155 187 89  ; ...
     128 100 162 ; ...
     75 172 198  ; ...
     247 150 70  ; ...
     206 185 102 ; ...
     156 176 132 ; ...
     107 177 201 ; ...
     101 133 207 ; ...
     126 107 201 ; ...
     163 121 187 ; ...
     83 84 138 ; ...
     67 128 134 ; ...
     160 77 163 ; ...
     196 101 45 ; ...
     139 93 61 ; ...
     92 146 181 ; ...
     110 160 176 ; ...
     204 175 10 ; ...
     141 137 164 ; ...
     116 133 96 ; ...
     158 146 115 ; ...
     126 132 141] / 255;  % color map
mymap(1:1:M,:) = A(mod([1:1:M]-1,24)+1,:) ;
% explicitly set the colormap's first entry to (white, pk 11/24/2021) black for background
mymap(1,:)=[1 1 1];
mymap(2,:)=[0 0 0];
% end Tommy, Mar 14, 2009
% -------------------

% % % % -------------------
% % % % begin Tommy, Apr 21, 2009
% % % % cells perimeter
% % % cellPerim = zeros(size(L,1),size(L,2)) ;
% % % for i1 = 1:1:max(L(:))
% % %    cellPerim = cellPerim + bwperim(L==i1) ;
% % % end
% % % % end Tommy, Apr 21, 2009
% % % % -------------------

% % get sequence of random integers in range [1,maxcolors-1]
% [s,I] = sort(rand(M-1,1));  
% % randomly reorder mymap color entries [2,maxcolors]
% mymap(2:end,:) = mymap(I+1,:);

% get sequence of random integers in range [1,maxcolors-1]
[s,I] = sort(rand(M-2,1));  
% randomly reorder mymap color entries [2,maxcolors]
mymap(3:end,:) = mymap(I+2,:);

if nargin>=2,
	rgb = 0.5 * ind2rgb(L2,mymap);
	bwscreen = double(bwscreen);
	bwscreen = 0.5 * bwscreen / maxmax(bwscreen);
	rgb(:,:,1) = rgb(:,:,1) + bwscreen;
	rgb(:,:,2) = rgb(:,:,2) + bwscreen;
	rgb(:,:,3) = rgb(:,:,3) + bwscreen;
	imshow(rgb); hold on ;
	outim = rgb;
else
	imshow(L2, mymap); hold on ;
    outim = L2;
end
% % % % -------------------
% % % % begin Tommy, Apr 21, 2009
% % % spy(cellPerim,'w'); hold off ; % cells perimeter
% % % % end Tommy, Apr 21, 2009
% % % % -------------------
    