%-----------------------------------------------------------%
%                                                           %
%   SEM Image Texture Analysis                              %
%   Eric Wilcox-Freeburg                                    %
%   Rev 0.2a    Last Modified 02/30/13                      %
%                                                           %
%   Notes: Histogram equalization should yield balanced     %
%   images.  Appears to be a gradient effect still across   %
%   the images, making the edges erroneously "smooth" due   %
%   to low contrast between values along those regions.     %
%   MATLAB help has an example of gradient reduction.       %
%                                                           %
%   Want to add pattern recognition to analyze sample for   %
%   patterns as opposed to box counting methods.            %
%                                                           %
%-----------------------------------------------------------%


clear all;
close all;

%% Point directory for batch processing
ftypes = {'*.jpg','JPEG';'*.bmp','Bitmap';'*.tif','TIF'};
[Filename, PATH, FIndx] = uigetfile(ftypes,'Select image:','C:\MATLAB');
RPF = [PATH,Filename];
switch FIndx;
    case 1;
        FTa = '*.jpg';
    case 2;
        FTa = '*.bmp';
    case 3;
        FTa = '*.tif';
end
cd(PATH);

items = dir(FTa);           % Index images of same file type
params = dir('*.txt');      % Identify parameter files if saved (JEOL format)


%%
OUTComp = [];                           % Preallocate compiled output matrix
for samp = 1:size(items);
    close all;
%% Load image (index by lit in cell array items().name)
    I = imread(items(samp).name);
    [D1, FileNm, D2] = fileparts(items(samp).name);
    
%% Histogram equalization
    if size(I,3)==3;
         I = .2989*I(:,:,1)+.5870*I(:,:,2)+.1140*I(:,:,3);
    end
    Ir = imresize(I, [NaN 2048]);
    Eim = histeq(Ir);
    BW1 = im2bw(Eim, 0.66);

    %%
    BWao = bwareaopen(BW1,8000);
    nhood = true(9);
    closeBWao = imclose(BWao,nhood);        % Smooth
    roughMask = imfill(closeBWao,'holes');  % Fill in holes in mask
    I2 = Ir;
    I2(~roughMask) = NaN;                   % Background subtracted raw image
    imshow(I2);

    %%
    button3 = questdlg(items(samp).name,'Skip?','Skip','Cont','Cont');  % Remove error volumes
%%
    switch button3
        case 'Skip'
            b3 = 1;
        case 'Cont'
            b3 = 0;
    end
    %%
    if b3 == 0;
        %% Remove extraneous features.    
            button2 = questdlg('Erroneous extra volume data?');  % Remove error volumes
            switch button2
                case 'Yes'
                    b2 = 1;
                case 'No'
                    b2 = 0;
                case 'Cancel'
                    b2 = 2;
            end

        %% Loop to create multiple ROIs for deletion
            delc = [];
            delr = [];
            while b2 == 1;
                imshow(roughMask);
                redc = [];
                redr = [];
                if b2 == 1;
                    erd = impoly;
                    pos = getPosition(erd);
                    c = pos(:,1);
                    r = pos(:,2);
                    BWi = (roipoly(roughMask,c,r));
                    [redc, redr] = find(BWi==1);
                end
                delc = [delc; redc];
                delr = [delr; redr];
                imshow(I2);
                    button2 = questdlg('Additional erroneous extra volume data?');  % Remove error volumes
                    switch button2
                        case 'Yes'
                            b2 = 1;
                        case 'No'
                            b2 = 0;
                        case 'Cancel'
                            b2 = 2;
                    end
            end

        %% subtract from roughMask
            zs = size(Ir);
            sr = size(delr);
            Mask1 = zeros(zs(1),zs(2));
            for erdc = 1:sr(1);
                roughMask(delc(erdc), delr(erdc)) = 0;  % Subtract exclusion area
            end
        %% construct new mask
            roughMask = bwareaopen(roughMask,8000); 
            I2(~roughMask) = NaN;                   % Background subtracted raw image
            imshow(I2);


        %%
        sigma1 = 10;
        sigma2 = 20;
        Igrad=localnormalize(I2,sigma1,sigma2);
        Imean=mean(Igrad(roughMask));
        imshow(Igrad)

        %% Analyze standard deviation of pixel-neighbors ~ texture
            S = stdfilt(Igrad, nhood);
            S = (S./Imean).*100;                    % RSD
            S2 = S>90;                            	% High RSD = high roughness?
            S3 = (~S2);                             % Inverts to give "smoothness" and removes otolith edge.
            S3(~roughMask) = 0;                     % Removes background
            imshow(S3);


        %% Calculate object statistics
            PxCount = sum(sum(S3));                 % # of "smooth" pixels
            STATS = regionprops(roughMask,'all');   % Calculates eccentricity
            AreaCount = sum(sum(roughMask));        % Area of ROI
            SmRat = PxCount/AreaCount;              % "Smoothness" ratio.
            Ecc = STATS(1).Eccentricity;            % Pull eccentricity data for fitted elipse.


        %% Save data to temporary matrix, append to compiled matrix
            OutTemp = {FileNm SmRat Ecc};
            OUTComp = [OUTComp; OutTemp];

        %% Save figures to files
            imwrite(I2,[FileNm '_bcksubt.bmp'],'bmp');  % Save background subtracted figure
            imwrite(S3,[FileNm '_proc.bmp'],'bmp');     % Save processed image
        %% Generate ellipse overlay
            EOv = fig('units','pixels','width',zs(2),'height',zs(1));
            imshow(S3)
            hold on
            phi = linspace(0,2*pi,50);
            cosphi = cos(phi);
            sinphi = sin(phi);

            for k = 1:length(STATS)
                xbar = STATS(k).Centroid(1);
                ybar = STATS(k).Centroid(2);

                a = STATS(k).MajorAxisLength/2;
                b = STATS(k).MinorAxisLength/2;

                theta = pi*STATS(k).Orientation/180;
                R = [ cos(theta)   sin(theta)
                    -sin(theta)   cos(theta)];

                xy = [a*cosphi; b*sinphi];
                xy = R*xy;

                x = xy(1,:) + xbar;
                y = xy(2,:) + ybar;

                plot(x,y,'r','LineWidth',2);
            end
            hold off
            %%
            export_fig(EOv, [FileNm '_overlay.bmp'], '-bmp');
    else
    end
end

%% Create header for output file
header = {'Sample' 'Smoothness Ratio' 'Eccentricity'};
OUTx = [header; OUTComp];
xlswrite('Summary.xls', OUTx, 1);