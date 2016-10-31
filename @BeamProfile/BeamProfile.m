classdef BeamProfile < Video
    %BeamProfile Class to handle beam profiles and their interaction with a spherical drop
    %   The class is supposed to help handling beam profiles as measured, for example, by a CCD beam
    %   profiler for pulsed lasers, i.e. measurements of the fluence (energy per area). It is a
    %   subclass of the Video class, extended by functions and properties required to work with beam
    %   profiles.
    %
    %   The functions for an actual beam profile are defined for a 2D matrix, either a dimensionless
    %   fluence of the order of one or a dimensional fluence in J/m^2. Use the idxFrames and
    %   idxChannels property to select certain frames and channels in the video that are then used
    %   to compute a single fluence matrix. That way, the raw data can be stored efficiently as
    %   video, for example as unit16 values. The video data is converted to a normalized fluence
    %   (fluenceNorm) taking into account a normalisation with a BASE value and a SCALE (stored in
    %   video2Norm, use estimateVideo2Norm to determine the value automatically),
    %              <normalized value> = (double(<video data>) - video2Norm(1)) / video2Norm(2).
    %   The dimensional fluence, for example in J/m^2, is optained by
    %             <dimensional value> = <normalized value> * norm2Abs,
    %   where norm2Abs can be set directly, otherwise it is computed to match a given mean energy
    %   per frame for the selected frames (just set the energy and norm2Abs is determined
    %   accordingly). Please note: norm2Abs is supposed to be the same for all frames in the video,
    %   i.e. the beam profiles need to be recorded with a constant relation between the pixel values
    %   and dimensional fluence, e.g. do not change the attenuation of the beam or gain of the CCD
    %   during the measurement. It is also assumed that all recording took place at a fixed position
    %   in space, i.e. the position of the CCD is the same for all frames.
    %
    %   Coordinate systems: the global video coordinate system defines a fixed XYZ system to
    %   describe the experimental arrangement. The positon of the CCD and its normal need to be
    %   given in this coordinate system, where z = 0 refers to the beam waist position of the laser.
    %   When talking about a beam profile or any image, x and y are the horizontal and vertical
    %   coordinate in the plane perpendicular to the laser beam path z (again, z = 0 is the position
    %   of the beam waist). The beam coordinate system xBeam and yBeam define a coordinate system,
    %   the origin of which is at the centroid of the beam calculated according to ISO11146-1. This
    %   coordinate system is valid for the current fluence (calculated based on selection by
    %   idxFrames and idxChannels).
    %
    %   Note: The class stores only properties to disk mandatory to recover the object, e.g. the ISO
    %   computation is repeated when an objected is loaded to disk since the results have not been
    %   stored to disk.
    %
    %-----------------------------------------------------------------------------------------------
    %   Copyright 2016 Alexander Ludwig Klein, alexludwigklein@gmail.com
    %
    %   Physics of Fluids, University of Twente
    %-----------------------------------------------------------------------------------------------
    
    %% Properties
    properties (GetAccess = public, SetAccess = public, Dependent = true)
        % idxFrames The index of the selected video frame(s) to use for the current beam profile (double)
        idxFrames
        % idxChannels The index of the selected video channel(s) to use for the current beam profile (double)
        idxChannels
        % energy The total energy in J of the current beam profile (double)
        energy
        % video2Norm Base value and scale to transform video data to normalized fluence (double)
        video2Norm
        % norm2Abs Calibration factor from normalized fluence to absolute fluence (double)
        norm2Abs
        % wavelength Wavelength of laser beam in m (double)
        wavelength
        % polarization Rotation of plane of polarization in rad relative to (yz), NaN for circular polarized light (double)
        %
        % Example: polarization = pi/4 means: the plane of polarization is rotated by pi/4 clockwise
        % compared to the horzontal image axis
        polarization
        % focalLength Effective focal length of the focussing optics in m (double)
        focalLength
        % zBeam Coordinate along z axis relative to beam centroid in Z, same as z and zCCD (double)
        zBeam
    end
    
    properties (GetAccess = public, SetAccess = protected, Dependent = true)
        % fluenceNorm Normalized fluence of current frames (matrix of double)
        fluenceNorm
        % fluence Absolute fluence (energy density in J/m^2) of current frames (matrix of double)
        fluence
        % xBeam Coordinates of pixels along beam x axis in m (nX float)
        xBeam
        % yBeam Coordinates of pixels along beam y axis in m (nY float)
        yBeam
        % xBeamGrid Grid of CCD for beam x in m (nX x nY float)
        xBeamGrid
        % yBeamGrid Grid of CCD for beam y in m (nX x nY float)
        yBeamGrid
        % zBeamGrid Grid of CCD for beam z in m (nX x nY float)
        zBeamGrid
        % ixMax X index of point of maximum fluence (double)
        ixMax
        % iyMax Y index of point of maximum fluence (double)
        iyMax
        % xMax Global x position in m of point of maximum fluence (double)
        xMax
        % yMax Global y position in m of point of maximum fluence (double)
        yMax
        % xBeamMax Beam x position in m of point of maximum fluence (float)
        xBeamMax
        % yBeamMax Beam y position in m of point of maximum fluence (float)
        yBeamMax
        % iso11146 Single beam profile properties according to ISO11146-1 (structure)
        iso11146
        % m1X First order moment in m in global x direction according to ISO11146-1, i.e. X centroid (float)
        m1X
        % m1Y First order moment in m in global y direction according to ISO11146-1, i.e. Y centroid (float)
        m1Y
        % m2X Second order moment in m^2 in x direction according to ISO11146-1 (float)
        m2X
        % m2Y Second order moment in m^2 in y direction according to ISO11146-1 (float)
        m2Y
        % m2XY Second order moment in m^2 in xy direction according to ISO11146-1 (float)
        m2XY
        % dX Beam width in m in x direction according to ISO11146-1 (float)
        dX
        % dY Beam width in m in y direction according to ISO11146-1 (float)
        dY
        % dR Beam diameter in m based on second order moments according to ISO11146-1 (float)
        dR
        % roit Rotation in rad of beam according to ISO11146-1 (float)
        rot
        % ell Ellipticity of beam according to ISO11146-1 (float)
        ell
        % maskBeam Mask for the beam based on ISO11146-1 properties (logical with size of fluenceNorm)
        maskBeam
        % maskInt Mask of interrogation area used in computations of ISO11146-1 properties (logical with size of fluenceNorm)
        maskInt
        % ecc Eccentricity of beam (float)
        ecc
        % dEff Effective beam diameter of a circular beam based on maskBeam (float)
        dEff
        % AFlatMax Effective area of a flat top with same peak fluence (float)
        AFlatMax
        % AFlatMean Effective area of a flat top with same mean fluence (float)
        AFlatMean
        % dFlatMax Effective diameter of a circular flat top with same peak fluence (float)
        dFlatMax
        % dFlatMean Effective diameter of a circular flat top with same mean fluence (float)
        dFlatMean
        % fluenceMax Maximum fluence within mask of beam in J/m^2 (float)
        fluenceMax
        % fluenceMean Mean fluence within mask of beam in J/m^2 (float)
        fluenceMean
        % fluenceMin Minimum fluence within mask of beam in J/m^2 (float)
        fluenceMin
    end
    
    properties (GetAccess = protected, SetAccess = protected, Transient = true)
        % p_iso11146 Storage for iso11146
        p_iso11146    = [];
        % p_energy Storage for energy
        p_energy      = [];
        % p_fluenceNorm Storage for fluenceNorm
        p_fluenceNorm = [];
        % p_xBeamGrid Storage for xBeamGrid
        p_xBeamGrid   = [];
        % p_yBeamGrid Storage for yBeamGrid
        p_yBeamGrid   = [];
        % p_zBeamGrid Storage for zBeamGrid
        p_zBeamGrid   = [];
        % p_ixMax Storage for ixMax
        p_ixMax       = [];
        % p_iyMax Storage for iyMax
        p_iyMax       = [];
        % p_fluenceMax Storage for fluenceMax
        p_fluenceMax  = [];
        % p_fluenceMean Storage for fluenceMean
        p_fluenceMean  = [];
        % p_fluenceMin Storage for fluenceMin
        p_fluenceMin  = [];
    end
    
    % the following properties are stored to disk and must be accessible for superclass
    properties (GetAccess = {?Video,?BeamProfile}, SetAccess = {?Video,?BeamProfile}, Transient = true)
        % p_video2Norm Storage for video2Norm
        p_video2Norm  = [];
        % p_norm2Abs Storage for norm2Abs
        p_norm2Abs    = [];
        % p_idxFrames Storage for idxFluence
        p_idxFrames   = 1;
        % p_idxChannels Storage for idxFluence
        p_idxChannels = 1;
        % p_wavelength Storage for wavelength
        p_wavelength  = NaN;
        % p_polarization Storage for polarization
        p_polarization= NaN;
        % p_focalLength Storage for focalLength
        p_focalLength = NaN;
    end
    
    %% Constructor, SET/GET
    methods
        function obj   = BeamProfile(filename,varargin)
            % BeamProfile Class constructor accepting one or more data/filename(s) and options
            %
            % The first input is the filename to one or multiple files as char or cellstr. The
            % options are parsed with MATLAB's input parser. The behaviour of the constructor is the
            % same as for a Video object with a bit of extension:
            %   * As subclass of the Video class it suports reading from TIF, AVI, MP4, MJ2 or its
            %     own and simple DAT format. All additional properties are stored in a MAT file of
            %     the same basename but extension '.mat'.
            %   * Compared to the video class the BeamProfiler class supports to read data from CSV
            %     files written by Thorlabs software for beam profiles. For each file one instance
            %     of the class is created, but memory mapping is not supported. Use the function
            %     convertCSV2VID to convert multiple CSV files to a single TIF, AVI, MP4 or MJ2.
            
            %
            % treat input of CSV files, set filename to empty for those elements and set their
            % actual values after the call to the class constructor in the supclass
            if nargin < 1, filename = []; end
            if ischar(filename), filename = {filename}; end
            if iscellstr(filename)
                bak              = filename;
                [~,~,ext]        = cellfun(@(x) fileparts(x),filename,'un',false);
                idxCSV           = ismember(ext,{'.csv','.CSV'});
                filename(idxCSV) = {''};
            else
                idxCSV = false;
            end
            %
            % redirect input to superclass constructor and add props to save to disk
            obj@Video(filename,varargin{:});
            %
            % read CSV files
            if any(idxCSV)
                idxCSV = find(idxCSV);
                for i = reshape(idxCSV,1,[])
                    dat                = obj.readCSV(bak{i});
                    [fp, fn]           = fileparts(bak{i});
                    obj(i).cdata.cdata = dat.data;          %#ok<AGROW>
                    setProperties(obj(i), dat);
                    obj(i).filename    = fullfile(fp,fn);   %#ok<AGROW>
                    obj(i).name        = 'CSV import';      %#ok<AGROW>
                end
            end
        end
        
        function         set.idxFrames(obj,value)
            if isa(value,'double') && min(round(value(:))) > 0 && max(round(value(:))) <= obj.nFrames
                value = reshape(round(value),1,[]);
                if ~isequal(value,obj.p_idxFrames)
                    resetUpdate(obj);
                    obj.p_idxFrames = value;
                end
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for idxFrames');
            end
        end
        
        function value = get.idxFrames(obj)
            if isempty(obj.p_idxFrames)
                obj.p_idxFrames = 1;
            end
            value = obj.p_idxFrames;
        end
        
        function         set.idxChannels(obj,value)
            if isa(value,'double') && min(round(value(:))) > 0 && max(round(value(:))) <= obj.nZ
                value = reshape(round(value),1,[]);
                if ~isequal(value,obj.p_idxChannels)
                    resetUpdate(obj);
                    obj.p_idxChannels = value;
                end
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for idxChannels');
            end
        end
        
        function value = get.idxChannels(obj)
            if isempty(obj.p_idxChannels)
                obj.p_idxChannels = 1;
            end
            value = obj.p_idxChannels;
        end
        
        function         set.energy(obj,value)
            if isa(value,'double') && isscalar(value)
                resetUpdate(obj);
                obj.p_norm2Abs = value/obj.pixresR.^2/(sum(obj.fluenceNorm(:)));
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for energy');
            end
        end
        
        function value = get.energy(obj)
            if isempty(obj.p_energy)
                obj.p_energy = obj.norm2Abs * sum(obj.fluenceNorm(:)) * obj.pixresR.^2;
            end
            value = obj.p_energy;
        end
        
        function         set.norm2Abs(obj,value)
            if isa(value,'double') && isscalar(value)
                resetUpdate(obj);
                obj.p_norm2Abs = value;
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for norm2Abs');
            end
        end
        
        function value = get.norm2Abs(obj)
            if isempty(obj.p_norm2Abs)
                if isempty(obj.p_energy), obj.p_energy = 1; end
                obj.p_norm2Abs = obj.p_energy/obj.pixresR.^2/(sum(obj.fluenceNorm(:)));
            end
            value = obj.p_norm2Abs;
        end
        
        function         set.video2Norm(obj,value)
            if isa(value,'double') && numel(value) == 2
                resetUpdate(obj);
                obj.p_video2Norm = reshape(value,1,2);
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for video2Norm');
            end
        end
        
        function value = get.video2Norm(obj)
            if isempty(obj.p_video2Norm)
                estimateVideo2Norm(obj);
            end
            value = obj.p_video2Norm;
        end
        
        function         set.wavelength(obj,value)
            if isa(value,'double') && isscalar(value)
                resetUpdate(obj);
                obj.p_wavelength = value;
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for wavelength');
            end
        end
        
        function value = get.wavelength(obj)
            value = obj.p_wavelength;
        end
        
        function         set.polarization(obj,value)
            if isa(value,'double') && isscalar(value)
                resetUpdate(obj);
                obj.p_polarization = value;
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for polarization');
            end
        end
        
        function value = get.polarization(obj)
            value = obj.p_polarization;
        end
        
        function         set.fluenceNorm(obj,value) %#ok<INUSD>
            error(sprintf('%s:Input',mfilename),['The normalized fluence cannot be set directly, ',...
                'change the cdata or the video2Norm property accordingly']);
        end
        
        function value = get.fluenceNorm(obj)
            if isempty(obj.p_fluenceNorm)
                % return the mean of selected frame(s) and channel(s)
                nSel  = numel(obj.idxFrames);
                value = mean(double(obj.cdata(:,:,obj.idxChannels,obj.idxFrames(1))),3)/nSel;
                for n = 2:nSel
                    value = value + mean(double(obj.cdata(:,:,obj.idxChannels,obj.idxFrames(n))),3)/nSel;
                end
                value = (value-obj.video2Norm(1))/obj.video2Norm(2);
                if obj.p_bufferData, obj.p_fluenceNorm = value; end
            else
                value = obj.p_fluenceNorm;
            end
        end
        
        function         set.fluence(obj,value) %#ok<INUSD>
            error(sprintf('%s:Input',mfilename),['The fluence cannot be set directly, ',...
                'change the cdata or the video2Norm property accordingly']);
        end
        
        function value = get.fluence(obj)
            value = obj.fluenceNorm * obj.norm2Abs;
        end
        
        function         set.focalLength(obj,value)
            if isa(value,'double') && isscalar(value) && (isnan(value) || value >= 0)
                resetUpdate(obj);
                obj.p_focalLength = value;
            else
                error(sprintf('%s:Input',mfilename),'Input not valid for focal length');
            end
        end
        
        function value = get.focalLength(obj)
            value = obj.p_focalLength;
        end
        
        function value = get.xBeam(obj)
            value = obj.x - obj.m1X;
        end
        
        function value = get.yBeam(obj)
            value = obj.y - obj.m1Y;
        end
        
        function value = get.zBeam(obj)
            value = obj.zCCD;
        end
        
        function         set.zBeam(obj,value)
            obj.zCCD = value;
        end
        
        function value = get.xBeamGrid(obj)
            if isempty(obj.p_xBeamGrid) || size(obj.p_xBeamGrid,1) ~= obj.nY || size(obj.p_xBeamGrid,2) ~= obj.nX
                [value, tmp] = meshgrid(obj.xBeam,obj.yBeam);
                if obj.p_bufferData
                    obj.p_xBeamGrid = value;
                    obj.p_yBeamGrid = tmp;
                end
            else
                value = obj.p_xBeamGrid;
            end
        end
        
        function value = get.yBeamGrid(obj)
            if isempty(obj.p_yBeamGrid) || size(obj.p_yBeamGrid,1) ~= obj.nY || size(obj.p_yBeamGrid,2) ~= obj.nX
                [tmp, value] = meshgrid(obj.xBeam,obj.yBeam);
                if obj.p_bufferData
                    obj.p_xBeamGrid = tmp;
                    obj.p_yBeamGrid = value;
                end
            else
                value = obj.p_yBeamGrid;
            end
        end
        
        function value = get.zBeamGrid(obj)
            if isempty(obj.p_zBeamGrid) || size(obj.p_zBeamGrid,1) ~= obj.nY || size(obj.p_zBeamGrid,2) ~= obj.nX
                value = repmat(cast(obj.zCCD,'like',obj.m1X),obj.nX,obj.nY);
                if obj.p_bufferData, obj.p_zBeamGrid = value; end
            else
                value = obj.p_zBeamGrid;
            end
        end
        
        function value = get.iyMax(obj)
            if isempty(obj.p_iyMax)
                % query obj.ixMax that sets obj.ixMax and obj.iyMax
                obj.ixMax;
            end
            value = obj.p_iyMax;
        end
        
        function value = get.ixMax(obj)
            if isempty(obj.p_ixMax)
                [~, ilin]                 = max(obj.fluenceNorm(:));
                [obj.p_iyMax,obj.p_ixMax] = ind2sub(size(obj.fluenceNorm),ilin);
            end
            value = obj.p_ixMax;
        end
        
        function value = get.xMax(obj)
            value = obj.x(obj.ixMax);
        end
        
        function value = get.yMax(obj)
            value = obj.y(obj.iyMax);
        end
        
        function value = get.xBeamMax(obj)
            value = obj.xBeam(obj.ixMax);
        end
        
        function value = get.yBeamMax(obj)
            value = obj.yBeam(obj.iyMax);
        end
        
        function value = get.iso11146(obj)
            if isempty(obj.p_iso11146)
                runISO11146(obj);
            end
            value = obj.p_iso11146;
        end
        
        function value = get.m1X(obj)
            value = obj.iso11146.m1X;
        end
        
        function value = get.m1Y(obj)
            value = obj.iso11146.m1Y;
        end
        
        function value = get.m2X(obj)
            value = obj.iso11146.m2X;
        end
        
        function value = get.m2Y(obj)
            value = obj.iso11146.m2Y;
        end
        
        function value = get.m2XY(obj)
            value = obj.iso11146.m2XY;
        end
        
        function value = get.dX(obj)
            value = obj.iso11146.dX;
        end
        
        function value = get.dY(obj)
            value = obj.iso11146.dY;
        end
        
        function value = get.dR(obj)
            value = obj.iso11146.dR;
        end
        
        function value = get.rot(obj)
            value = obj.iso11146.rot;
        end
        
        function value = get.ell(obj)
            value = obj.iso11146.ell;
        end
        
        function value = get.maskBeam(obj)
            value = obj.iso11146.maskBeam;
        end
        
        function value = get.maskInt(obj)
            value = obj.iso11146.maskInt;
        end
        
        function value = get.ecc(obj)
            value = obj.iso11146.ecc;
        end
        
        function value = get.dEff(obj)
            value = obj.iso11146.dEff;
        end
        
        function value = get.AFlatMax(obj)
            value = obj.iso11146.AFlatMax;
        end
        
        function value = get.AFlatMean(obj)
            value = obj.iso11146.AFlatMean;
        end
        
        function value = get.dFlatMax(obj)
            value = sqrt(4/pi*obj.iso11146.AFlatMax);
        end
        
        function value = get.dFlatMean(obj)
            value = sqrt(4/pi*obj.iso11146.AFlatMean);
        end
        
        function value = get.fluenceMax(obj)
            if isempty(obj.p_fluenceMax)
                obj.p_fluenceMax = max(obj.fluence(obj.maskBeam(:)));
            end
            value = obj.p_fluenceMax;
        end
        
        function value = get.fluenceMean(obj)
            if isempty(obj.p_fluenceMean)
                obj.p_fluenceMean = mean(obj.fluence(obj.maskBeam(:)));
            end
            value = obj.p_fluenceMean;
        end
        
        function value = get.fluenceMin(obj)
            if isempty(obj.p_fluenceMin)
                obj.p_fluenceMin = min(obj.fluence(obj.maskBeam(:)));
            end
            value = obj.p_fluenceMin;
        end
    end
    
    %% Methods for various class related tasks
    methods (Access = public, Hidden = false)
        function value      = isdefault(obj)
            %isdefault Tests if object is the default object, i.e. contains no data worth storing,
            % please note: this is not a 100% fail-safe test
            
            value = isdefault@Video(obj);
            for i = reshape(find(value),1,[])
                value(i) = (isnan(obj(i).wavelength) && isnan(obj(i).polarization) && ...
                    isnan(obj(i).focalLength) &&...
                    isscalar(obj(i).idxFrames) && obj(i).idxFrames == 1 && ...
                    isscalar(obj(i).idxChannels) && obj(i).idxChannels == 1);
            end
        end
        
        function              runISO11146(obj,varargin)
            % runISO Runs the ISO11146-1 calculation and sets property accordingly
            
            %
            % process one-by-one
            if numel(obj) > 1
                for i = 1:numel(obj)
                    runISO11146(obj(i),varargin{:});
                end
                return;
            end
            %
            % check input
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            %  initial threshold in realitve fluence to find beam
            opt.addParameter('thres', 0.1, ...
                @(x) isfloat(x) && isscalar(x));
            opt.parse(varargin{:});
            opt = opt.Results;
            %
            % find initial interrogation area (rectangle which holds everything above a certain
            % threshold), rectangle is increased as suggested in general by ISO.
            int2beam = 2;   % size of interrogation area compared to beam size
            nPix     = round(0.05*sqrt(obj.nX^2+obj.nX^2));
            tmp      = obj.fluenceNorm ;
            tmp      = tmp > opt.thres*max(tmp(:));
            tmp      = bwareaopen(tmp,round(nPix));
            tmp      = imdilate(tmp,strel('disk',nPix));
            stats    = regionprops(tmp,'SubarrayIdx','Area');
            if numel(stats) ~= 1
                warning(sprintf('%s:ISO11146',mfilename),['Fluence data seems to contain %d ',...
                    'beams and, therfore, is not unambiguous. The largest object is used to ',...
                    'find the beam. Please check!'], numel(stats));
                if numel(stats) < 1
                    stats = struct;
                    stats.SubarrayIdx{2} = 1:obj.nX;
                    stats.SubarrayIdx{1} = 1:obj.nY;
                else
                    [~,idxObj] = max([stats.Area]);
                    stats      = stats(idxObj);
                end
            end
            idxX     = stats.SubarrayIdx{2};
            idxY     = stats.SubarrayIdx{1};
            [xGridObj, yGridObj] = meshgrid(obj.x,obj.y); % local computation for faster access
            F        = obj.fluenceNorm(idxY,idxX);
            xg       = xGridObj(idxY,idxX);
            yg       = yGridObj(idxY,idxX);
            % perform iterative process to determine beam propertie
            doAgain = true;
            counter = 2;
            nMax    = 100;
            res     = Inf;
            while doAgain
                % 1st order moments (eq. 1 & eq. 2 in ISO)
                E   = sum(F(:));
                m1X = sum(F(:).*xg(:))./E;                          %#ok<PROPLC>
                m1Y = sum(F(:).*yg(:))./E;                          %#ok<PROPLC>
                % 2nd order moments (eq. 3 & eq. 4 in ISO)
                m2X  = sum(F(:).*(xg(:)- m1X).^2)./E;               %#ok<PROPLC>
                m2Y  = sum(F(:).*(yg(:)- m1Y).^2)./E;               %#ok<PROPLC>
                m2XY = sum(F(:).*(xg(:)- m1X).*(yg(:)- m1Y))./E;    %#ok<PROPLC>
                % beam diameter (eq. 15 to 23 in ISO)
                if abs(m2X-m2Y)/(abs(m2X)+abs(m2Y)) > 1e-2          %#ok<PROPLC>
                    rot   = 0.5 * atan(2*m2XY/(m2X-m2Y));           %#ok<PROPLC>
                    gamma = sign(m2X-m2Y);                          %#ok<PROPLC>
                    tmp1  = m2X+m2Y;                                %#ok<PROPLC>
                    tmp2  = gamma.*((m2X-m2Y).^2+4*m2XY.^2).^0.5;   %#ok<PROPLC>
                else
                    rot  = sign(m2XY)*pi/4;                         %#ok<PROPLC>
                    tmp1 = m2X+m2Y;                                 %#ok<PROPLC>
                    tmp2 = 2*abs(m2XY);                             %#ok<PROPLC>
                end
                dX  = 8^0.5 .*(tmp1+tmp2).^0.5;                     %#ok<PROPLC>
                dY  = 8^0.5 .*(tmp1-tmp2).^0.5;                     %#ok<PROPLC>
                dR  = 8^0.5 .*(tmp1).^0.5;                          %#ok<PROPLC>
                % ellipticity
                ell        = min(dX,dY)./max(dX,dY);                %#ok<PROPLC>
                % check if ellipticity is converged
                if counter >= nMax
                    warning(sprintf('%s:ISO11146',mfilename),['Calculation of parameters is not ',...
                        'converged after %d iterations. Please check!'], nMax);
                    doAgain = false;
                elseif abs(ell-res) < 0.005 %#ok<PROPLC>
                    doAgain = false;
                else
                    % create new interrogation zone
                    idxX    = find(obj.x >= m1X - int2beam/2 * dX & obj.x <= m1X + int2beam/2 * dX); %#ok<PROPLC>
                    idxY    = find(obj.y >= m1Y - int2beam/2 * dY & obj.y <= m1Y + int2beam/2 * dY); %#ok<PROPLC>
                    F       = obj.fluenceNorm(idxY,idxX);
                    xg      = xGridObj(idxY,idxX);
                    yg      = yGridObj(idxY,idxX);
                    res     = ell; %#ok<PROPLC>
                    counter = counter + 1;
                end
            end
            % add results to property
            obj.p_iso11146.m1X     = m1X;   %#ok<PROPLC>
            obj.p_iso11146.m1Y     = m1Y;   %#ok<PROPLC>
            obj.p_iso11146.m2X     = m2X;   %#ok<PROPLC>
            obj.p_iso11146.m2Y     = m2Y;   %#ok<PROPLC>
            obj.p_iso11146.m2XY    = m2XY;  %#ok<PROPLC>
            obj.p_iso11146.dX      = dX;    %#ok<PROPLC>
            obj.p_iso11146.dY      = dY;    %#ok<PROPLC>
            obj.p_iso11146.dR      = dR;    %#ok<PROPLC>
            obj.p_iso11146.rot     = rot;   %#ok<PROPLC>
            obj.p_iso11146.ell     = ell;   %#ok<PROPLC>
            % eccentricity
            obj.p_iso11146.ecc     = sqrt(max(dX,dY).^2-min(dX,dY).^2)./max(dX,dY);                                             %#ok<PROPLC>
            %
            % add a few properties not directly mentioned in ISO11146, e.g. masks to easily find
            % the beam: make mask for beam and interrogation zone, mask for beam is an ellipse.
            maskInt            = false(size(obj.fluenceNorm));                                                                  %#ok<PROPLC>
            maskInt(idxY,idxX) = true;                                                                                          %#ok<PROPLC>
            maskBeam           = ((xGridObj-m1X).*cos(rot)+(yGridObj-m1Y).*sin(rot)).^2./(dX/2).^2 + ...                        %#ok<PROPLC>
                ((xGridObj-m1X).*sin(rot)-(yGridObj-m1Y).*cos(rot)).^2./(dY/2).^2 <= 1;                                         %#ok<PROPLC>
            obj.p_iso11146.maskInt   = maskInt;                                                                                 %#ok<PROPLC>
            obj.p_iso11146.maskBeam  = maskBeam;                                                                                %#ok<PROPLC>
            obj.p_iso11146.dEff      = cast(sqrt(sum(maskBeam(:))*obj.pixresR^2*4/pi),'like',m1X);                              %#ok<PROPLC>
            obj.p_iso11146.AFlatMax  = sum(sum(obj.fluenceNorm(maskBeam)))/max(max(obj.fluenceNorm(maskBeam)))*obj.pixresR^2;   %#ok<PROPLC>
            obj.p_iso11146.AFlatMean = sum(sum(obj.fluenceNorm(maskBeam)))/mean(mean(obj.fluenceNorm(maskBeam)))*obj.pixresR^2; %#ok<PROPLC>
            % check baseline in interrogation zone without beam
            if isnan(dR)                                                                                                        %#ok<PROPLC>
                warning(sprintf('%s:ISO11146',mfilename),'Computation for ISO11146 leads to NaN');
            else
                maskNoise    = maskInt & ~maskBeam;                                                                             %#ok<PROPLC>
                noise2signal = sum(sum(obj.fluenceNorm(maskNoise)))/sum(sum(obj.fluenceNorm(maskBeam)));                        %#ok<PROPLC>
                if noise2signal > 0.20
                    warning(sprintf('%s:ISO11146',mfilename),['Noise to signal ratio seems ',...
                        'to be too high: %.2e. Please check!'],noise2signal);
                end
            end
        end
        
        function              estimateVideo2Norm(obj,varargin)
            % estimateVideo2Norm Estimates and sets video2Norm property to normalise video data
            % based on the currently selected frames
            %
            % The algorithm to characterize a laser beam is based on second order moments that are
            % strongly influenced by a false baseline in the fluence. Therefore, a background level
            % is estimated here. It is important to keep negative noise, since it leads to noise
            % cancellation when the scond order moments are computed. Finally the data is scaled to
            % be of the order of one.
            
            %
            % process one-by-one
            if numel(obj) > 1
                for i = 1:numel(obj)
                    estimateVideo2Norm(obj(i),varargin{:});
                end
                return;
            end
            if isnan(obj), obj.video2Norm = NaN(1,2); return; end
            %
            % check input
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % initial threshold in realitve fluence to find beam
            opt.addParameter('thres', 0.1, ...
                @(x) isfloat(x) && isscalar(x));
            % true false whether to show the beam found in the image data
            opt.addParameter('plot', false, ...
                @(x) islogical(x) && isscalar(x));
            opt.parse(varargin{:});
            opt = opt.Results;
            %
            % get fluence without any normalisation
            bak               = obj.p_video2Norm;
            obj.p_video2Norm  = [0 1];
            obj.p_fluenceNorm = [];
            value             = obj.fluenceNorm;
            obj.p_video2Norm  = bak;
            obj.p_fluenceNorm = [];
            %
            % try to find a beam
            v2n      = [1 1];
            nPix     = round(0.05*sqrt(size(value,1)^2+size(value,2)^2));
            tmp      = value > opt.thres*max(value(:));
            tmp      = bwareaopen(tmp,round(nPix));
            tmp      = imdilate(tmp,strel('disk',nPix));
            stats    = regionprops(tmp,'SubarrayIdx');
            if numel(stats) ~= 1
                warning(sprintf('%s:Input',mfilename),['Data matrix seems to contain %d beams and, ',...
                    'therfore, is not unambiguous, please check!'],numel(stats));
            end
            if sum(tmp(:)) > 0.9*size(value,1)*size(value,2)
                warning(sprintf('%s:Input',mfilename),['Beam seems to cover almost complete data matrix (>90%%) and ',...
                    'baseline may not be correct, please check!']);
            end
            
            %
            % estimate noise level and scale
            v2n(1)       = sum(value(~tmp(:)))/sum(~tmp(:));
            noise2signal = v2n(1)/(sum(value(tmp(:)))/sum(tmp(:)));
            if noise2signal > 0.1
                warning(sprintf('%s:Input',mfilename),['Noise to signal ratio seems to be too ',...
                    'high: %.2e, please check!'],noise2signal);
            end
            v2n(2)           = max(value(:) - v2n(1));
            resetUpdate(obj);
            obj.p_video2Norm = v2n;
            %
            % plot beam(s) in fluence
            if opt.plot
                img = value - min(value(:));
                img = img./max(img(:));
                img = Video.imgColorize(img,tmp);
                fig = figure('Name','Noise estimation'); %#ok<NASGU>
                imshow(img);
                title(sprintf('Normalized image and %d beam(s), video2Norm = [%s]',...
                    numel(stats),num2str(obj.video2Norm)));
            end
        end
        
        function obj        = sort(obj)
            % sort Sorts array of objects according to the position along the propagation axis z
            
            [~,idx] = sort([obj.zCCD]);
            obj     = obj(idx);
        end
        
        function out        = average(obj,varargin)
            % average Averages array of beam profiles and returns new object with one average, the
            % cdata class of the new object is double.
            
            %
            % parse and check input
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % quantity to average
            opt.addParameter('quantity', 'fluence', ...
                @(x) ischar(x) && ismember(x,{'fluence', 'fluenceNorm'}));
            opt.parse(varargin{:});
            opt = opt.Results;
            %
            % check number of pixel, pixel resolution, etc.
            checkPropertyEqual(obj,{'pixresX' 'pixresY' 'wavelength' 'polarization' 'xCCD' 'yCCD' 'focalLength'},1);
            checkPropertyEqual(obj,{'nX' 'nY' 'norm'},2);
            tmp = [obj(:).energy]; tmp = (max(tmp) - min(tmp))/max(tmp);
            if tmp > 0.05
                warning(sprintf('%s:Input',mfilename),['Averaging is performed for beam profiles ',...
                    'with a variation in energy of %.2g %%'],tmp*100);
            end
            %
            % make a new object
            out              = BeamProfile;
            setProperties(out,obj);
            out.name         = sprintf('Average of %d beam profile(s)',numel(obj));
            out.device       = sprintf('Average of %d beam profile(s)',numel(obj));
            out.comment      = sprintf('Result of averaging %d beam profile(s)',numel(obj));
            out.userdata     = [];
            out.date         = datetime('now');
            data             = cat(3,obj(:).(opt.quantity));
            data             = data./max(data(:));
            data             = mean(data,3);
            data             = data/max(data(:));
            out.cdata.cdata  = data;
            out.idxFrames    = 1;
            out.idxChannels  = 1;
            out.video2Norm   = [0 1];
            out.energy       = mean([obj(:).energy]);
        end
        
        function out        = flyBeam(obj,varargin)
            % flyBeam Creates video object with the current fluence of each object as single frame,
            % the output can be used to fly along the beam by 'out.play'
            
            %
            % parse and check input
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % quantity to fly through
            opt.addParameter('quantity', 'fluence', ...
                @(x) ischar(x) && ismember(x,{'fluence', 'fluenceNorm'}));
            opt.parse(varargin{:});
            opt = opt.Results;
            %
            % check number of pixel, pixel resolution, etc.
            checkPropertyEqual(obj,{'pixresX' 'pixresY' 'wavelength' 'polarization' 'xCCD' 'yCCD'},1);
            checkPropertyEqual(obj,{'nX' 'nY'},2);
            tmp = [obj(:).energy]; tmp = (max(tmp) - min(tmp))/max(tmp);
            if tmp > 0.05
                warning(sprintf('%s:Input',mfilename),['Concatenation is performed for beam ',...
                    'profiles with a variation in energy of %.2g %%'],tmp*100);
            end
            %
            % base output on first object
            out             = obj(1).copy;
            out.device      = sprintf('Concatenation of %d beam profiles',numel(obj));
            out.comment     = {sprintf('Result of combining %d beam profiles',numel(obj))};
            out.userdata    = [];
            out.date        = datetime('now');
            data            = cat(3,obj(:).(opt.quantity));
            data            = data./max(data(:));
            out.cdata.cdata = permute(data,[1 2 4 3]);
            out.energy      = mean([obj(:).energy]);
        end
        
        function varargout  = plotProfile(obj,varargin)
            % plot Plots the laser beam profile(s) and returns handles to graphic objects in output structure
            %
            % Options in varargin are passed to an input parser, see code
            %
            
            nargoutchk(0,1);
            %
            % use input parser to process options
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % show warning if there are too many objects to plot
            opt.addParameter('maxFigure', 10, ...
                @(x) isnumeric(x) && isscalar(x));
            % property to show: fluence or fluenceNorm
            opt.addParameter('showProp', 'fluence', ...
                @(x) ischar(x) && ismember(x,{'fluence','fluenceNorm'}));
            % property to show: fluence or fluenceNorm
            opt.addParameter('showInfo', false, ...
                @(x) islogical(x) && isscalar(x));
            % number of points to plot some objects, e.g. ellipse
            opt.addParameter('nPlot', 100, ...
                @(x) isnumeric(x) && isscalar(x));
            % true/false whether to reuse existing figure with correct tag when one object is shown
            opt.addParameter('plotReuse', true, ...
                @(x) islogical(x) && isscalar(x));
            opt.parse(varargin{:});
            opt = opt.Results;
            nColumn = 2+opt.showInfo;
            %
            % start plotting profiles one by one
            if numel(obj) > 1
                if numel(obj) > opt.maxFigure
                    button = questdlg(sprintf('This will plot %d beam profiles in new figures. Continue?',numel(obj)), ...
                        sprintf('%s:plotProfile',mfilename), 'Yes','No','No');
                    if strcmp(button,'No'); return; end
                end
                opt.plotReuse = false;
                for i = 1:numel(obj)
                    out(i) = plotProfile(obj(i), opt); %#ok<AGROW>
                end
                if nargout > 0, varargout = {out}; end
                return;
            elseif numel(obj) < 1
                error(sprintf('%s:Input',mfilename),'At least one beam profiles are required to plot')
            end
            %
            % plot single profile and prepare figure
            if opt.showInfo
                pos = [50, 50, 1400, 750];
            else
                pos = [50, 50, 1000, 750];
            end
            out.fig = [];
            if opt.plotReuse
                fig = findall(groot,'type','figure', 'Tag','BeamProfile_plotProfile');
                if ~isempty(fig) && isvalid(fig(1))
                    out.fig = fig(1);
                    delete(out.fig.Children);
                end
            end
            if isempty(out.fig)
                out.fig = figure('Name',sprintf('Beam profile from %s', datestr(obj.date)),...
                    'Position',pos,'Visible','off','HandleVisibility','callback',...
                    'Tag','BeamProfile_plotProfile');
            end
            out.ax = gobjects(5,1);
            %
            % prepare some data for plotting
            switch opt.showProp
                case 'fluence'
                    strTitle = sprintf('Fluence (J/cm^2): min = %.2g, mean = %.2g, max = %.2g',...
                        1e-4*min(obj.(opt.showProp)(:)), 1e-4*mean(obj.(opt.showProp)(:)), 1e-4*max(obj.(opt.showProp)(:))); %#ok<NASGU>
                    str2 = 'fluence (J/cm^2)';
                    data = 1e-4*obj.fluence;
                case 'fluenceNorm'
                    strTitle = sprintf('Normalized fluence: min = %.2g, mean = %.2g, max = %.2g',...
                        min(obj.(opt.showProp)(:)), mean(obj.(opt.showProp)(:)), max(obj.(opt.showProp)(:))); %#ok<NASGU>
                    str2 = 'normalized fluence';
                    data = obj.fluenceNorm;
            end
            strX = ['profile x, global ' obj.xStr ' (mm)']; strY =['profile y, global ' obj.yStr ' (mm)'];
            %
            % image plot
            out.ax(1) = subplot(2,nColumn,1,'NextPlot','Add','Parent',out.fig);
            out.image = imagesc(obj.x([1 end])*1e3,obj.y([1 end])*1e3, data,'Parent',out.ax(1));
            colorbar(out.ax(1),'northoutside');
            if obj.xDir > 0, set(out.ax(1),'XDir','normal');
            else,            set(out.ax(1),'XDir','reverse');
            end
            if obj.yDir > 0, set(out.ax(1),'YDir','reverse');
            else,            set(out.ax(1),'YDir','normal');
            end
            plot(out.ax(1), [min(obj.x) max(obj.x)]*1e3, [obj.y(obj.iyMax) obj.y(obj.iyMax)]*1e3,...
                'Color','black','Linestyle','--','DisplayName','Max');
            h = plot(out.ax(1), [obj.x(obj.ixMax) obj.x(obj.ixMax)]*1e3, [min(obj.y) max(obj.y)]*1e3,...
                'Color','black','Linestyle','--','DisplayName','Max');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            xlabel(out.ax(1),strX); ylabel(out.ax(1),strY); % title(out.ax(1),strTitle);
            % add ellipse of beam based on parametric equation and principal axis
            xE = @(t,scale) obj.m1X + scale*obj.dX/2 * cos(obj.rot) * cos(t) - scale*obj.dY/2 * sin(obj.rot) * sin(t);
            yE = @(t,scale) obj.m1Y + scale*obj.dX/2 * sin(obj.rot) * cos(t) + scale*obj.dY/2 * cos(obj.rot) * sin(t);
            t  = linspace(0,2*pi,opt.nPlot);
            plot(out.ax(1), xE(t,1)*1e3, yE(t,1)*1e3, 'Color','red','Linestyle','-',...
                'DisplayName','Beam');
            h = plot(out.ax(1), xE([0 pi],1.5)*1e3, yE([0 pi],1.5)*1e3, 'Color','red','Linestyle','-',...
                'DisplayName','Beam');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h = plot(out.ax(1), xE([pi/2 1.5*pi],1.5)*1e3, yE([pi/2 1.5*pi],1.5)*1e3, ...
                'Color','red','Linestyle','-','DisplayName','Beam');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            %
            % line plots at maxima
            out.ax(2) = subplot(2,nColumn,2,'NextPlot','Add','YDir','reverse','Parent',out.fig);
            plot(out.ax(2), data(:,obj.ixMax), obj.y*1e3, 'Color','black','DisplayName','Fluence at max');
            ylabel(out.ax(2), strY); xlabel(out.ax(2), str2);
            out.ax(3) = subplot(2,nColumn,3+opt.showInfo,'NextPlot','Add','XDir','reverse','Parent',out.fig);
            plot(out.ax(3), obj.x*1e3, data(obj.iyMax,:), 'Color','black','DisplayName','Fluence at max');
            xlabel(out.ax(3), strX); ylabel(out.ax(3), str2);
            % indicate position of maximum
            plot(out.ax(2), [0 data(obj.iyMax,obj.ixMax)], [obj.y(obj.iyMax) obj.y(obj.iyMax)]*1e3, ...
                'Color','black','Linestyle','--','DisplayName','Max');
            plot(out.ax(3), [obj.x(obj.ixMax) obj.x(obj.ixMax)]*1e3, [0 data(obj.iyMax,obj.ixMax)],...
                'Color','black','Linestyle','--','DisplayName','Max');
            % show mean value on perimeter of beam
            mask  = bwperim(obj.maskBeam);
            meanF = mean(mean(data(mask)));
            plot(out.ax(2), [meanF meanF], [min(obj.y) max(obj.y)]*1e3, ...
                'Color','red','DisplayName','Fluence at beam edge');
            plot(out.ax(3), [min(obj.x) max(obj.x)]*1e3, [meanF meanF], ...
                'Color','red','DisplayName','Fluence at beam edge');
            hLinkAx(1) = linkprop(out.ax([1 2]),'ylim');
            hLinkAx(2) = linkprop(out.ax([1 3]),'xlim');
            setappdata(out.fig,'hLinkAx',hLinkAx);
            %
            % 3D plot, might make panning crash in image plot
            out.ax(4) = subplot(2,nColumn,4+opt.showInfo,'Parent',out.fig);
            surf(out.ax(4),obj.xGrid*1e3, obj.yGrid*1e3, data,'LineStyle','none','FaceColor','Texture','EdgeColor','none');
            axis(out.ax(4), 'vis3d');
            if obj.xDir > 0, set(out.ax([3 4]),'XDir','normal');
            else,            set(out.ax([3 4]),'XDir','reverse');
            end
            if obj.yDir > 0, set(out.ax([2 4]),'YDir','reverse');
            else,            set(out.ax([2 4]),'YDir','normal');
            end
            xlabel(out.ax(4), strX);
            ylabel(out.ax(4), strY);
            zlabel(out.ax(4), str2);
            % information on beam
            if opt.showInfo
                out.ax(5) = subplot(2,nColumn,[3;6],'Parent',out.fig);
                axis(out.ax(5),'off');
                hInfo = text(0,1.1,obj.info,'Parent',out.ax(5),'HorizontalAlignment','left',...
                    'VerticalAlignment','top', 'FontName','FixedWidth','Interpreter','none');
            else
                out.ax(5) = [];
            end
            %
            % final thoughts
            set(out.ax(1),'XTickLabel', []);
            set(out.ax(1),'XLabel', []);
            set(out.ax(2),'YTickLabel', []);
            set(out.ax(2),'YLabel', []);
            set(out.ax,   'LineWidth',2,'Fontsize',14);
            set(findall(out.fig,'type','text'),'fontSize',14);
            set(out.ax(1),'XLim', 1e3*[min(obj.x) max(obj.x)]);
            set(out.ax(1),'YLim', 1e3*[min(obj.y) max(obj.y)]);
            if opt.showInfo
                pos1 = [0.05 0.48 0.33 0.45];
            else
                pos1 = [0.1 0.48 0.6 0.45];
            end
            set(out.ax(1),'Position', pos1);
            pos1 = get(out.ax(1),'Position');
            pos2 = [pos1(1)+pos1(3)+0.02 pos1(2) 0.2 pos1(4)];
            pos3 = [pos1(1) 0.05 pos1(3) pos1(2)-0.05-0.04];
            pos4 = [pos1(1)+pos1(3)+0.05 0.05 pos2(3)-0.05 pos3(4)-0.05];
            pos5 = [pos4(1)+pos4(3)+0.03 0 1-pos4(1)-pos4(3) pos2(2)+pos2(4)];
            set(out.ax(2),'Position', pos2);
            set(out.ax(3),'Position', pos3);
            set(out.ax(4),'Position', pos4);
            if opt.showInfo
                set(hInfo,'fontSize',12,'FontName','FixedWidth','Interpreter','none');
                set(out.ax(5),'Position', pos5);
            end
            axis(out.ax(1),'equal');
            set(out.fig,'Visible','on');
            if nargout > 0
                varargout = {out};
            end
        end
        
        function varargout  = placeDrop(obj,varargin)
            %placeDrop Puts a drop in the beam and estimates the fluence on its surface and the
            % amount transmitted into the surface taking into account the cosine effect and Fresnel
            % equations. The function also returns also intermediate results to have everything in
            % one structure for later usage. Optional input to this function is read by inputparser,
            % have a look at the code for available options
            %
            % Idea of computation: The fluence of the laser beam, given by the beam profile, is
            % measured at the drop position. Each pixel of this planar beam profile is interpreted
            % as a ray and for each ray it is possible to determine the position and incident angle
            % where the ray hits the drop in case the drop is present in the laser beam path.
            % Assuming the optical axis of the beam is aligned with the drop coordinate system,
            % this raytracing is a matter of matrix operations and can be performed quickly.
            %
            % The output structure:
            %       beamprofile: The beam profile used for the computation, which may be the result
            %                    of an interpolation or average (BeamProfile)
            %
            %   The following fields of the structure describe the geometry of the problem in a drop
            %   coordinate system (drop's center-of-mass is at x=y=0, laser beam waist is at z=0):
            %                 n: Refractive index of the surrounding and drop, respectively (double)
            %                R0: Initial drop radius in m (double)
            %               xy0: Drop xy position relative to beam (double)
            %              rMax: Maximum radius on CCD from which rays can hit the drop, i.e. an
            %                    aperture that indicates which rays do not reach the drop (double)
            %            phiMax: Maximum phi angle that can be reached by a ray (double)
            %                 x: x coordinate of each ray where it hits the drop (double)
            %                 y: y coordinate of each ray where it hits the drop (double)
            %                 r: r coordinate of each ray where it hits the drop (double)
            %              mask: Mask of rays hitting the drop in the image of the fluence, please
            %                    note: the mask may be larger or smaller than the actual drop in
            %                    case the rays are not assumed to be parallel (logical)
            %               phi: Angle in spherical coordinate system where each ray hits the drop,
            %                    runs from 0 (center of the drop where laser hits dead on) to pi/2
            %                    (outer edge of the drop when viewed in laser direction z) to pi
            %                    (back of the drop) (double)
            %             theta: Angle in spherical coordinate system where each ray hits the drop,
            %                    runs from -pi to pi (double)
            %             phiIn: Incident angle at the drop surface, 0 means perpendicular impact
            %                    (center of drop) and pi/2 means parallel, in case of a parallel
            %                    beam phiIn equals phi (double)
            %            phiFix: Correction added to phi to account for incoming angle of rays due
            %                    to focussing when calculating phiIn, i.e. phiIn = phi + phiFix (double)
            %
            %   The following fields describe the input fluence and several convolution matrices to
            %   account for the cosine effect and the fresnel equations:
            %       fluenceNorm: The normalized, planar fluence at the drop's center-of-mass
            %                    position as measured by the beam profiler, in case the option
            %                    'perfectBeam' is set to true, this is the fluence of a perfect flat
            %                    top beam profile of equal energy and size (double)
            %          norm2Abs: A scale to transform the normalized fluence to a dimensional value in J/m^2 (double)
            %           convCOS: Effect of the drop's curvature on the fluence, i.e. the 'cosine
            %                    effect' called for the earth and the solar energy reaching the
            %                    curved surface of the earth, runs from 0 to 1 (double)
            %       convFresnel: Effect of Fresnel equations on the fluence, runs from 0 to 1 (double)
            %            convCF: The transmitted normalized fluence into the drop, i.e:
            %                       convCF = fluenceNorm .* convCOS .* convFresnel (double)
            %         convCFAbs: The transmitted absolute fluence into the drop, i.e:
            %                    convCFAbs = fluenceNorm .* convCOS .* convFresnel *  norm2Abs (double)
            %
            %   An axi-symmetric representation is computed and fitted to a gaussian and a gaussian
            %   convoluted with a cosine yielding an analytical expression for the energy
            %   transmitted into the drop as function of phi:
            %          parGauss: Parameters for an axi-symmetric description of the transmitted
            %                    fluence into the drop at its surface based on a gaussian relation:
            %                             T =  par(1) * exp(-phi.^2/(2*par(2)^2));
            %       parCOSGauss: Parameters for an axi-symmetric description of the transmitted
            %                    fluence into the drop at its surface based on a cos * gaussian
            %                    relation:
            %                             T = par(1) * cos(phiIn) .* exp(-phi.^2/(2*par(2)^2));
            %                    In this relation, phiIn is the incident angle, which is
            %                           phiIn = phi, for a parallel beam
            %                           phiIn = phi + atan(sin(phi)./(z0/R0-cos(phi)). for z0 < 0
            %                           phiIn = phi + atan(sin(phi)./(z0/R0+cos(phi)). for z0 > 0
            %                    The value z0/R0 is stored as third element in the output
            %            linPhi: All phi of the drop pixels
            %             linCF: All convCF of the drop pixels
            %            fitPhi: Consolidated linPhi used for the fit (otherwise biased weights in fit)
            %             fitCF: Consolidated linCF used for the fit (otherwise biased weights in fit)
            %
            %   Based on the matrices, several energies and fluences can be calculated:
            %          e2D_beam: Total energy of the beam in J (for the sake of completeness) (double)
            %        e2D_planar: Total energy in projected area in front of the drop in J (double)
            %         e2D_trans: Total energy transmitted into the drop in J (double)
            % e2D_relative_beam: Fraction of the total beam energy transmitted to the drop (double)
            % e2D_relative_drop: Fraction of the beam energy in the projected area transmitted to
            %                    into the drop (e.g. 1 means that all energy in pi R0^2  was
            %                    transmitted to the drop, which is the case when Fresnel equations
            %                    are neglected) (double)
            %    f2D_planar_max: Maximum planar (not wrapped around the drop) fluence in J/m^2 (double)
            %   f2D_planar_mean: Mean planar fluence in J/m^2 (double)
            %    f2D_planar_min: Minimum planar fluence in J/m^2 (double)
            %      f2D_curv_max: Maximum fluence on drop surface (wrapped around the drop) in J/m^2 (double)
            %     f2D_curv_mean: Mean fluence on drop surface in J/m^2 (double)
            %      f2D_curv_min: Minimum fluence on drop surface in J/m^2 (double)
            %     f2D_trans_max: Maximum fluence transmitted through drop surface (wrapped around the drop) in J/m^2 (double)
            %    f2D_trans_mean: Mean fluence transmitted through drop surface in J/m^2 (double)
            %     f2D_trans_min: Minimum fluence transmitted through drop surface in J/m^2 (double)
            %    f2D_e2D_planar: e2D_planar devided by pi R0^2 in J/m^2, same as f2D_planar_mean (double)
            %     f2D_e2D_trans: e2D_trans devided by pi R0^2 in J/m^2, different from f2D_trans_mean (double)
            %      f2D_e2D_beam: Mean planar fluence of the complete beam based on dR in ISO11146-1 (for the sake of completeness) (double)
            %
            %   The figure(s) and axes in case the result are plotted:
            %               fig: The figure object(s) when plotting (graphics, figure)
            %                ax: The axes object(s) when plotting (graphics, figure)
            
            nargoutchk(0,1);
            %
            % use input parser to process options
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % drop radius in m
            opt.addParameter('R0', 1e-3, ...
                @(x) isnumeric(x) && isscalar(x) && x > 0);
            % refractive index for the surrounding of the drop and the liquid of the drop
            opt.addParameter('n', [1 1.33], ...
                @(x) isnumeric(x) && isvector(x) && numel(x) == 2);
            % position of drop along z axis, takes the mean of given profiles if it is empty
            % otherwise a beam profile is looked for measured at this position, an interpolation is
            % performed if no beam profile is found
            opt.addParameter('z', [], ...
                @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
            % position of drop in fluence plane (XY plane) relative to beam center
            opt.addParameter('xy', [0 0], ...
                @(x) isnumeric(x) && isvector(x) && numel(x) == 2);
            % true/false whether to just fit to theoretical prediction of a perfect flat top beam
            % profile, i.e. the actual fluence of the beam is assumed to be the same everywhere.
            % This can be used to illustrate the effect of fresnel equation, for example. The
            % fluence is matched to get the same energy as the actual beam (assuming the perfect
            % beam fills the same mask)
            opt.addParameter('perfectBeam', false, ...
                @(x) islogical(x) && isscalar(x));
            % true/false whether to perform the calculation assuming a paralle beam, this neglects
            % the position of the drop relative to the beam waist, same as assumuning the beam waist
            % is at infinity
            opt.addParameter('parallelBeam', false, ...
                @(x) islogical(x) && isscalar(x));
            % imresize fluence by certain factor for better resolution
            opt.addParameter('resize', 1, ...
                @(x) isnumeric(x) && isscalar(x) && x > 0);
            % true/false whether to plot results with default settings in placeDropPlot
            opt.addParameter('plot', true, ...
                @(x) islogical(x) && isscalar(x));
            opt.parse(varargin{:});
            opt = opt.Results;
            %
            % find beam profile to work with
            if isempty(opt.z), opt.z = mean([obj.z]); end
            idx = find(abs([obj.z]-opt.z) < 1e-6,1,'first');
            if isempty(idx)
                if numel(obj) < 2 || ~(opt.z > min([obj.z]) && opt.z < max([obj.z]))
                    error(sprintf('%s:Input',mfilename),['Interpolation of beam profiles can not be ',...
                        'performed for given drop position, either not enough profiles or position ',...
                        'is out of bounds, please check!']);
                end
                bp = interpolate(obj, opt.z);
            else
                % Option 1: copy complete beam profile
                % bp = obj(idx).copy;
                % % copy results of iso calculation to speed up
                % bp.p_iso11146 = obj(idx).p_iso11146;
                %
                % Option 2: Make a shallow copy and set a new fluence fluence (requires less memory)
                tmpE           = obj(idx).energy;
                flu            = obj(idx).fluenceNorm;
                bp             = obj(idx).copyShallow;
                bp.cdata.cdata = flu;
                bp.idxFrames   = 1;
                bp.idxChannels = 1;
                bp.video2Norm  = [0 1];
                bp.name        = ['Copy of beam profile: ' obj(idx).name];
                bp.device      = ['Copy of beam profile: ' obj(idx).device];
                bp.comment     = cat(1,{'Copy of beam profile: '},reshape(obj(idx).comment,[],1));
                bp.date        = datetime('now');
                bp.userdata    = [];
                bp.reset;
                % make that energy is correct
                bp.energy = tmpE;
            end
            %
            % check if model is valid
            if ~opt.parallelBeam && abs(opt.R0/bp.z) > 1
                error(sprintf('%s:Input',mfilename),['Distance of the drop to the laser beam ',...
                    'waist is of the order of the initial drop radius, which may lead to errors ',...
                    'in the calculation. Use the parallel beam assumption, implement a model for ',...
                    'gaussian beam propagation and check the Rayleigh length or ... please check!']);
            end
            if sqrt(sum(opt.xy.^2))/bp.dR > 0.25
                warning(sprintf('%s:Input',mfilename),['Offset of the drop from the laser beam axis ',...
                    'is of the order of the beam radius, which leads to errors in the ray tracing. ',...
                    'Implement an asymmetric model for the ray tracing or ... please check!']);
            end
            if bp.z > 0
                warning(sprintf('%s:Input',mfilename),['Drop is placed after the beam waist, i.e. ',...
                    'z > 0, which is rather unlikely for a high enery beam ... please check!']);
            end
            %
            % * make a perfect flat top the same size as the actual beam
            % * determine the normalized input fluence to match the total energy in the actual beam
            % * make one pixel value in the perfect beam slightly greater to match the maximum in
            %   the actual beam and ease the ISO11146-1 computation.
            if opt.perfectBeam
                tmpE           = bp.energy;
                flu            = single(bp.maskBeam);
                [~,idxX]       = min(abs(bp.x-bp.m1X));
                [~,idxY]       = min(abs(bp.y-bp.m1Y));
                flu(idxY,idxX) = 1+1e-6 .* max(flu(:));
                bp.cdata.cdata = flu;
                bp.video2Norm  = [0 1];
                bp.energy      = tmpE;
                bp.reset;
            end
            if opt.resize ~= 1
                tmpE           = bp.energy;
                pixres         = bp.pixres;
                % Option 1: make a copy of selected frames and resize
                % flu            = bp.cdata.cdata(:,:,:,obj.idxFrames);
                % flu            = imresize(flu,opt.resize);
                % bp.cdata.cdata = flu;
                % bp.idxFrames   = 1:size(flu,4);
                %
                % Option 2: make a copy of just the fluence and resize (requires less memory)
                flu            = bp.fluenceNorm;
                flu            = imresize(flu,opt.resize);
                bp.cdata.cdata = flu;
                bp.idxFrames   = 1;
                bp.idxChannels = 1;
                bp.video2Norm  = [0 1];
                bp.pixres      = pixres/opt.resize;
                bp.reset;
                % make that energy is the same
                bp.energy      = tmpE;
            end
            out.beamprofile = bp;
            fluenceNorm     = bp.fluenceNorm; %#ok<PROPLC>
            %
            % * compute x and y, i.e. the coordinates on the CCD in a cartesian drop coordinate system
            % * compute phi and theta, i.e the coordinates on the CCD in a spherical drop coordinate
            %   system, phi is the same as: acos(sqrt(opt.R0^2-x.^2-y.^2)/opt.R0), which should come
            %   from the definition of spherical coordinates assuming r = R0 (for the points on the
            %   drop surface)
            % * compute the mask of the drop in the fluence image, where the mask includes all rays
            %   that reach the drop surface
            % * compute incident angle phiIn on drop surface (0 means perpendicular and pi/2 is
            %   parallel to the drop surface), instead of linear optics one could estimate the
            %   curvature of the wavefront based on a gaussian beam.
            %
            % posDrop is the distance of the drop to the beam waist, positive when drop is placed
            % after the beam waist, so z = 0 is the position of the beam waist, start computing x0,
            % etc., i.e. coordinates assuming a parallel beam
            posDrop  = bp.z;
            x0       = bp.xBeamGrid - opt.xy(1);
            y0       = bp.yBeamGrid - opt.xy(2);
            r0       = sqrt(x0.^2+y0.^2);
            mask0    = r0 <= opt.R0;
            % theta is not affected by the correction in r since we assume the drop's center-of-mass
            % is on the optical axis of the laser beam
            theta = atan2(y0,x0);
            if opt.parallelBeam
                % drop in parallel beam, illuminated for phi <= pi/2
                x      = x0;
                y      = y0;
                r      = r0;
                rMax   = opt.R0;
                phiMax = pi/2;
                mask   = mask0;
                phi    = asin(r0/opt.R0);
                phiFix = zeros(size(x));
            else
                % correct radius for a non-parallel beam
                r2z = (r0/posDrop).^2;
                if posDrop < 0
                    % drop before beam waist, illuminated for phi > pi/2
                    r1          = r0 .* (1./(1+r2z) + sqrt(1./(1+r2z).^2 - (1-r2z.*(opt.R0./r0).^2)./(1+r2z)));
                    phi         = NaN(size(mask0));
                    phi(mask0)  = asin(r1(mask0)/opt.R0);       % rays hitting the front of the drop
                    phi(~mask0) = pi-asin(r1(~mask0)/opt.R0);   % rays hitting the back of the drop
                    % compute the maximum radius of a ray that hits the drop
                    r0Max       = abs(tan(asin(abs(opt.R0/posDrop))) * posDrop);
                    mask        = r0<=r0Max;
                    phiMax      = pi/2+asin(abs(opt.R0/posDrop));
                else
                    % drop after beam waist, illuminated for phi < pi/2
                    r1    = r0 .* (1./(1+r2z) - sqrt(1./(1+r2z).^2 - (1-r2z.*(opt.R0./r0).^2)./(1+r2z)));
                    phi   = asin(r1/opt.R0);
                    r0Max = abs(tan(asin(abs(opt.R0/posDrop))) * posDrop);
                    mask  = r0<=r0Max;
                    phiMax= acos(abs(opt.R0/posDrop));
                end
                % compute correction for incident angle
                phiFix = atan(r0/posDrop); % atan(sin(phi)./(posDrop/opt.R0 - cos(phi)));
                % re-compute x and y
                x    = cos(theta) .* r1;
                y    = sin(theta) .* r1;
                r    = r1;
                rMax = r0Max;
            end
            if ~any(mask(:))
                error(sprintf('%s:Input',mfilename),'Beam does not seem to hit drop at all, please check');
            end
            % set rays not hitting the drop to default values
            phi(~mask)   = 0;
            x(~mask)     = 0;
            y(~mask)     = 0;
            r(~mask)     = 0;
            phiIn        = phi + phiFix;
            phiIn(~mask) = 0;
            % store current results in output structure
            out.n           = opt.n;
            out.R0          = opt.R0;
            out.xy0         = opt.xy;
            out.rMax        = rMax;
            out.phiMax      = phiMax;
            out.x           = x;
            out.y           = y;
            out.r           = r;
            out.mask     	= mask;
            out.phi         = phi;
            out.theta       = theta;
            out.phiIn       = phiIn;
            out.phiFix      = phiFix;
            out.fluenceNorm = fluenceNorm; %#ok<PROPLC>
            out.norm2Abs    = bp.norm2Abs;
            % * compute effect on fluence at the drop surface due to the cosine effect, as it is
            %   called for the earth and the solar energy reaching the curved surface of the earth.
            %   Here, we also assume a parallel incident radiation, i.e. the paraxial approximation.
            %   Note: this could be changed by estimating the curvature of the wavefront based on a
            %   gaussian beam.
            convCOS        = cos(phiIn);
            convCOS(~mask) = 0;
            out.convCOS    = convCOS;
            % * compute effect on fluence at the drop surface due to Fresnel equations
            Rs = @(n1,n2,w) abs((n1.*cos(w)-n2.*sqrt(1-(n1./n2.*sin(w)).^2))./(n1.*cos(w)+n2.*sqrt(1-(n1./n2.*sin(w)).^2))).^2;
            Rp = @(n1,n2,w) abs((n1.*sqrt(1-(n1./n2.*sin(w)).^2)-n2.*cos(w))./(n1.*sqrt(1-(n1./n2.*sin(w)).^2)+n2.*cos(w))).^2;
            Ts = @(n1,n2,w) 1-Rs(n1,n2,w);
            Tp = @(n1,n2,w) 1-Rp(n1,n2,w);
            R  = @(n1,n2,w) (Rs(n1,n2,w) + Rp(n1,n2,w))/2;
            T  = @(n1,n2,w) 1-R(n1,n2,w);
            if ~isnan(bp.polarization)
                % decompose field in S and P, fluence ~ field^2
                convFresnel = Ts(opt.n(1),opt.n(2),phiIn) .* sin(theta+bp.polarization).^2 + ...
                    Tp(opt.n(1),opt.n(2),phiIn) .* cos(theta+bp.polarization).^2;
            else
                % assuming circular polarized light
                convFresnel = T(opt.n(1),opt.n(2),phiIn);
            end
            convFresnel(~mask) = 0;
            out.convFresnel    = convFresnel;
            % * combine effects
            convCF        = out.fluenceNorm.*convCOS.*convFresnel;
            out.convCF    = convCF;
            out.convCFAbs = out.norm2Abs*convCF;
            %
            % compute a axi-symmetric representation and fit a prediction, either a gaussian or a
            % gaussion convoluted by a cosine
            linPhi       = reshape(phi(mask),1,[]);
            linCF        = reshape(convCF(mask),1,[]);
            [linPhi,idx] = sort(linPhi);
            linCF        = linCF(idx);
            % consolidate data to a mean curve for the fit
            [fitPhi, fitCF] = consolidator(linPhi,linCF,@mean,(pi/2)/100);
            % fit gaussian
            par0         = [max(linCF),     pi/8];
            parLB        = [0, 0];
            parUB        = [1.5*max(linCF), pi];
            options      = optimoptions('lsqcurvefit');
            options.Display = 'off';
            parGauss     = lsqcurvefit(@BeamProfile.myGaussian,par0,fitPhi,fitCF,parLB,parUB,options);
            out.parGauss = parGauss;
            % fit cos * gaussian
            par0         = [max(linCF),     pi/8];
            parLB        = [0, 0];
            parUB        = [1.5*max(linCF), pi*1e3];
            options      = optimoptions('lsqcurvefit');
            options.Display = 'off';
            if opt.parallelBeam
                parCOSGauss = lsqcurvefit(@BeamProfile.myCOSGaussian,par0,fitPhi,fitCF,parLB,parUB,options);
            else
                parCOSGauss = lsqcurvefit(@(par,xdata) BeamProfile.myCOSGaussianNonParallel([par bp.z/opt.R0],xdata),...
                    par0,fitPhi,fitCF,parLB,parUB,options);
                parCOSGauss = [parCOSGauss bp.z/opt.R0];
            end
            out.parCOSGauss = parCOSGauss;
            % add data used for the fit to output
            out.linPhi   = linPhi;
            out.linCF    = linCF;
            out.fitPhi   = fitPhi;
            out.fitCF    = fitCF;
            %
            % compute how much energy is given to the drop, note: the calculated normalized fluence
            % (convCF) is calculated on the drop surface. Therefore, the area of the pixels must
            % be altered (they are larger due to the curvature). This means to invers the effect of
            % convCOS, or just leave it out here. The only effect considered here that can reduce
            % the energy absorbed by the drop is Fresnel equation, the cosine effect changes the
            % fluence on the drop surface but not the total energy given to the drop.
            convE = convFresnel.*fluenceNorm; %#ok<PROPLC>
            if abs(bp.energy - sum(fluenceNorm(:)) * bp.norm2Abs * bp.pixresR.^2) > eps %#ok<PROPLC>
                warning(sprintf('%s:Input',mfilename),['Something went wrong in the energy computation, ',...
                    'energies do not match, please check!']);
            end
            out.e2D_beam          = bp.energy;
            out.e2D_planar        = sum(fluenceNorm(mask(:))) * bp.norm2Abs * bp.pixresR.^2; %#ok<PROPLC>
            out.e2D_trans         = sum(convE(:)) * bp.norm2Abs * bp.pixresR.^2;
            out.e2D_relative_beam = out.e2D_trans / out.e2D_beam;
            out.e2D_relative_drop = out.e2D_trans / out.e2D_planar;
            % compute fluence values
            out.f2D_planar_max   = max(fluenceNorm(mask(:)))   * bp.norm2Abs; %#ok<PROPLC>
            out.f2D_planar_mean  = mean(fluenceNorm(mask(:)))  * bp.norm2Abs; %#ok<PROPLC>
            out.f2D_planar_min   = min(fluenceNorm(mask(:)))   * bp.norm2Abs; %#ok<PROPLC>
            out.f2D_curv_max     = max(convCOS(mask(:))  .* fluenceNorm(mask(:))) * bp.norm2Abs; %#ok<PROPLC>
            out.f2D_curv_mean    = mean(convCOS(mask(:)) .* fluenceNorm(mask(:))) * bp.norm2Abs; %#ok<PROPLC>
            out.f2D_curv_min     = min(convCOS(mask(:))  .* fluenceNorm(mask(:))) * bp.norm2Abs; %#ok<PROPLC>
            out.f2D_trans_max    = max(convCF(mask(:)))   * bp.norm2Abs;
            out.f2D_trans_mean   = mean(convCF(mask(:)))  * bp.norm2Abs;
            out.f2D_trans_min    = min(convCF(mask(:)))   * bp.norm2Abs;
            out.f2D_e2D_planar   = out.e2D_planar/(pi * opt.R0^2);
            out.f2D_e2D_trans    = out.e2D_trans/(pi * opt.R0^2);
            out.f2D_e2D_beam     = bp.energy/(pi/4*bp.dR^2);
            %
            % plot and return results
            if opt.plot, [out.fig, out.ax] = BeamProfile.placeDropPlot(out); end
            if nargout > 0, varargout = {out}; end
        end
        
        function [out, F]   = interpolate(obj, newPos, varargin)
            % interpolate Interpolates array of beam profiles at new z position(s) and returns new
            % objects and the interpolant for the fluence. The interpolant uses pixel coordinates
            % for the x and y direction and the dimensional values for the z coordinate. It also
            % sorts the object along the z axis. The cdata class of the new object(s) is double;
            
            %
            % parse and check input
            nObj = numel(obj);
            if  nObj < 2
                error(sprintf('%s:Input',mfilename),'Interpolation requires at least two beam profiles for linear interpolation.');
            end
            if ~(isnumeric(newPos) && isvector(newPos))
                error(sprintf('%s:Input',mfilename),'Please specify the position of the interpolation as vector with z values.');
            end
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % method of interpolation
            opt.addParameter('method', 'linear', ...
                @(x) ischar(x) && ismember(x,{'linear','nearest','next','previous','pchip','cubic','spine'}));
            % method of extrapolation
            opt.addParameter('extrapolationMethod', 'linear', ...
                @(x) ischar(x) && ismember(x,{'linear','nearest','next','previous','pchip','cubic','spine'}));
            % quantity to interpolate
            opt.addParameter('quantity', 'fluence', ...
                @(x) ischar(x) && ismember(x,{'fluence', 'fluenceNorm'}));
            opt.parse(varargin{:});
            opt = opt.Results;
            %
            % check number of pixel, pixel resolution, etc.
            checkPropertyEqual(obj,{'pixresX' 'pixresY' 'wavelength' 'polarization' 'xCCD' 'yCCD' 'focalLength'},1);
            checkPropertyEqual(obj,{'nX' 'nY' 'norm'},2);
            tmp = [obj(:).energy]; tmp = (max(tmp) - min(tmp))/max(tmp);
            if tmp > 0.05
                warning(sprintf('%s:Input',mfilename),['Interpolation is performed for beam profiles ',...
                    'with a variation in energy of %.2g %%'],tmp*100);
            end
            %
            % sort object
            obj = sort(obj);
            %
            % create a 3D interpolant
            tmp = [[obj(:).xCCD]  [obj(:).yCCD]];
            tmp = sum((tmp - repmat(mean(tmp,1),size(tmp,1),1)).^2,2)^0.5;
            if tmp > eps
                % XY plane may differ among beam profiles: create a full grid, this options is
                % disabled for now and an error is thrown. It needs more checking, e.g. to flip the
                % grid for strictly monotonic increasing values, etc.
                % F = griddedInterpolant(cat(obj(:).xGrid,3),cat(obj(:).yGrid,3),cat(obj(:).zGrid,3),flu,...
                %     opt.method, opt.extrapolationMethod);
                error(sprintf('%s:Input',mfilename),['Arbitrary grids are currently not well supported ',...
                    'and the options is disabled, please extend the code :)']);
            else
                % XY plane is the same for all beam profiles, as should be the x and y vectors: use
                % compact representation for a XY grid in pixel coordinates and z in dimensional
                % form as given by the object(s). The fluence is normalized here to the global
                % maximum of all beam profiles, the rescale is done when each new object's energy is
                % set.
                flu      = cat(3,obj(:).(opt.quantity));
                flu      = flu./max(flu(:));
                mypixres = mean(cat(1,obj(:).pixres),1);
                myx      = (1:size(flu,2))*mypixres(1);
                myy      = (1:size(flu,1))*mypixres(2);
                F        = griddedInterpolant({myy(:),myx(:),[obj(:).z]},flu,...
                    opt.method, opt.extrapolationMethod);
            end
            %
            % create new beam profile(s) and interpolate
            nNew   = numel(newPos);
            newPos = {myy, myx, newPos(:)};
            out    = BeamProfile(nNew);
            for i = 1:nNew
                tmp    = newPos;
                tmp{3} = tmp{3}(i);
                setProperties(out(i),obj);
                out(i).cdata.cdata  = F(tmp);       % fluence could be rescaled to be within 0 to 1
                out(i).idxFrames    = 1;
                out(i).idxChannels  = 1;
                out(i).video2Norm   = [0 1];
                out(i).z            = tmp{3}(i);
                out(i).device       = sprintf('Interpolation - %s',opt.method);
                out(i).comment      = sprintf('Interpolated beam profile by %s method',opt.method);
                out(i).name         = sprintf('Interpolated beam profile by %s method',opt.method);
                out(i).energy       = mean([obj(:).energy]);
                out(i).userdata     = [];
                out(i).date         = datetime('now');
                out(i).pixres       = mypixres;
                out(i).norm         = obj(1).norm;
            end
        end
        
        function varargout  = plotPropagation(obj,varargin)
            % plotPropagation Fits a hyperbolic function to several beam profiles that are recorded
            % along z axis as described in ISO11146-1. Based on the fit, several beam parameters can
            % be determined, such as beam waist, etc. The results are added as userdata to first
            % object for later use and, finally, plotted.
            %
            % Options in varargin are passed to an input parser, see code
            %
            % Result is a structure with beam parameters based on dX, dY or dR:
            % 2Fit Fit of squared diameter
            %   z0 Location of beam waist
            %   d0 Diameter of beam waist
            %  div Beam divergence
            %   zR Rayleigh length of waist
            %   M2 M^2
            %
            
            nargoutchk(0,1);
            %
            % use input parser to process options
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % true/false whether to plot result
            opt.addParameter('plot', true, ...
                @(x) islogical(x) && isscalar(x));
            % number of points to plot some objects, e.g. ellipse
            opt.addParameter('nPlot', 100, ...
                @(x) isnumeric(x) && isscalar(x));
            opt.parse(varargin{:});
            opt = opt.Results;
            % check input
            if numel(obj) < 3
                error(sprintf('%s:Input',mfilename),'At least three beam profiles are required to perform a fit at all')
            elseif numel(obj) < 10
                warning(sprintf('%s:Input',mfilename),'At least 10 beam profiles are required to perform a good fit according to ISO11146-1')
            end
            %
            % check wavelength and polarization
            checkPropertyEqual(obj,{'wavelength' 'polarization'},1);
            % check stability of beam position
            tmp = [[obj(:).m1X]' [obj(:).m1Y]'];
            tmp = max(sqrt(sum((tmp - repmat(mean(tmp,1),size(tmp,1),1)).^2, 2)))/mean([obj.dR]);
            if tmp > 0.1
                warning(sprintf('%s:Input',mfilename),['Beam position seems to be unstable (maximum deviation from ',...
                    'average beam position relative to average beam size: %e %%'],tmp*100)
            end
            % check energy stability of beam
            tmp = [obj(:).energy]; tmp = (max(tmp)-min(tmp))/mean(tmp);
            if tmp > 0.05
                warning(sprintf('%s:Input',mfilename),['Beam energy seems to be unstable ',...
                    '(shot-to-shot variation: %.2g %%)'],tmp*100)
            end
            % fit each beam diameter squared to polynominal of degree 2
            res.zCCD  = [obj(:).zCCD];
            addFit('dX');
            addFit('dY');
            addFit('dR');
            addFit('dEff');
            addFit('dFlatMax');
            addFit('dFlatMean');
            % fit fluence
            res.fluenceMin      = [obj(:).fluenceMin];
            res.fluenceMin_Fit  = polyfit(res.zCCD,1./[obj(:).fluenceMin],2);
            res.fluenceMean      = [obj(:).fluenceMean];
            res.fluenceMean_Fit = polyfit(res.zCCD,1./[obj(:).fluenceMean],2);
            res.fluenceMax      = [obj(:).fluenceMax];
            res.fluenceMax_Fit  = polyfit(res.zCCD,1./[obj(:).fluenceMax],2);
            % warn when overwriting unknown data
            if ~(isempty(obj(1).userdata) || (isstruct(obj(1).userdata) && isequal(fieldnames(obj(1).userdata),fieldnames(res))))
                warning(sprintf('%s:Input',mfilename),'Replacing unknown userdata of first object with results of this function');
            end
            obj(1).userdata = res;
            %
            % Plot results
            if ~opt.plot
                return;
            end
            %
            % plot single profile and prepare figure
            out.fig = figure('Name',sprintf('Beam propagation based on %d beam profiles', numel(obj)),...
                'Position',[50, 50, 1300, 800],'Visible','off','HandleVisibility','callback');
            out.ax = gobjects(5,1);
            % plot diameters
            out.ax(1) = subplot(3,4,1:3,'NextPlot','Add','Parent',out.fig);
            ylabel(out.ax(1), 'diamter (mm)');
            xlabel(out.ax(1), 'laser beam axis z (mm)');
            mycolor   = distinguishable_colors(6);
            z         = linspace(min(res.zCCD)-0.1*(max(res.zCCD)-min(res.zCCD)),max(res.zCCD)+0.1*(max(res.zCCD)-min(res.zCCD)),opt.nPlot);
            i         = 1;
            plotFitD(out.ax(1),'dX');
            plotFitD(out.ax(1),'dY');
            plotFitD(out.ax(1),'dR');
            plotFitD(out.ax(1),'dEff');
            plotFitD(out.ax(1),'dFlatMax');
            plotFitD(out.ax(1),'dFlatMean');
            legend(out.ax(1),'show');
            plotWaist(out.ax(1));
            ylim = get(out.ax(1),'YLim');
            set(out.ax(1), 'YLim', [0 ylim(2)]);
            % plot relative change of energy, rotation and properties of ellipse
            out.ax(2) = subplot(3,4,5:7,'NextPlot','Add','Parent',out.fig);
            ylabel(out.ax(2), 'check stability');
            xlabel(out.ax(2), 'laser beam axis z (mm)');
            mycolor   = distinguishable_colors(4);
            str = 'energy'; tmp = [obj(:).(str)]; i = 1;
            plot(out.ax(2), 1e3*res.zCCD, tmp/mean(tmp), 'Displayname', str, ...
                'Linewidth', 2, 'Linestyle','none', 'Marker', 'o', 'Color', mycolor(i,:)); i = i + 1;
            str = 'rot'; tmp = [obj(:).(str)];
            plot(out.ax(2), 1e3*res.zCCD, tmp/pi, 'Displayname', str, ...
                'Linewidth', 2, 'Linestyle','none', 'Marker', 'o', 'Color', mycolor(i,:)); i = i + 1;
            str = 'ell'; tmp = [obj(:).(str)];
            plot(out.ax(2), 1e3*res.zCCD, tmp, 'Displayname', str, ...
                'Linewidth', 2, 'Linestyle','none', 'Marker', 'o', 'Color', mycolor(i,:)); i = i + 1;
            str = 'ecc'; tmp = [obj(:).(str)];
            plot(out.ax(2), 1e3*res.zCCD, tmp, 'Displayname', str, ...
                'Linewidth', 2, 'Linestyle','none', 'Marker', 'o', 'Color', mycolor(i,:)); i = i + 1;
            legend(out.ax(2),'show');
            % ylim = get(out.ax(2),'YLim');
            % set(out.ax(2), 'YLim', [0 ylim(2)]);
            % plot fluence
            out.ax(3) = subplot(3,4,9:11,'NextPlot','Add','Parent',out.fig);
            ylabel(out.ax(3), 'fluence (J/cm^2)');
            xlabel(out.ax(3), 'laser beam axis z (mm)');
            mycolor   = distinguishable_colors(3);
            i = 1;
            plotFitF(out.ax(3), 'fluenceMin');
            plotFitF(out.ax(3), 'fluenceMean');
            plotFitF(out.ax(3), 'fluenceMax');
            legend(out.ax(3),'show');
            % plot variation in position of beam in xy plane
            out.ax(4) = subplot(3,4,4,'NextPlot','Add','Parent',out.fig);
            xlabel(out.ax(4), 'beam stability y/dY');
            ylabel(out.ax(4), 'beam stability x/dX');
            set(out.ax(4),'YDir','reverse');
            set(out.ax(4),'XDir','reverse');
            plot(out.ax(4), ([obj(:).m1Y]-mean([obj(:).m1Y]))./[obj(:).dY], ([obj(:).m1X]-mean([obj(:).m1X]))./[obj(:).dX], 'Displayname', 'exp', ...
                'Linewidth', 2, 'Linestyle','none', 'Marker', 'o', 'Color', 'black');
            xlim = get(out.ax(4),'XLim');
            if max(abs(xlim)) < 1
                set(out.ax(4),'XLim',[-1 1]);
            end
            ylim = get(out.ax(4),'YLim');
            if max(abs(ylim)) < 1
                set(out.ax(4),'YLim',[-1 1]);
            end
            % plot info
            out.ax(5) = subplot(3,4,8,'Parent',out.fig);
            strInfo = {''};
            strInfo{2} = sprintf('|      Beam propagation      |');
            strInfo{1} = repmat('-',1,numel(strInfo{2}));
            strInfo{3} = strInfo{1};
            fn   = fieldnames(res);
            idx1 = find(cellfun(@(x) numel(x) > 2 && strcmp(x(end-2:end),'_M2'),fn));
            idx2 = find(cellfun(@(x) numel(x) > 2 && strcmp(x(end-3:end),'_div'),fn));
            idx3 = find(cellfun(@(x) numel(x) > 2 && strcmp(x(end-2:end),'_d0'),fn));
            idx4 = find(cellfun(@(x) numel(x) > 2 && strcmp(x(end-2:end),'_z0'),fn));
            strInfo{end+1} = sprintf('             M^2 | DIV');
            for k = 1:numel(idx1)
                if strcmp(fn{idx1(k)}(1:end-3),fn{idx2(k)}(1:end-4))
                    strInfo{end+1} = sprintf([' %9s: %4.2g | %4.2e',char(176)],fn{idx1(k)}(1:end-3),res.(fn{idx1(k)}),180/pi*res.(fn{idx2(k)})); %#ok<AGROW>
                else
                    error(sprintf('%s:Input',mfilename),'This should not happen');
                end
            end
            strInfo{end+1} = [' ' repmat('-',1,numel(strInfo{2})-2) ' '];
            strInfo{end+1} = sprintf('              d0');
            for k = 1:numel(idx3)
                strInfo{end+1} = sprintf(' %9s: %4.2e mm ',fn{idx3(k)}(1:end-3),1e3*res.(fn{idx3(k)})); %#ok<AGROW>
            end
            strInfo{end+1} = [' ' repmat('-',1,numel(strInfo{2})-2) ' '];
            strInfo{end+1} = sprintf('              z0');
            for k = 1:numel(idx4)
                strInfo{end+1} = sprintf(' %9s: %4.2e mm ',fn{idx4(k)}(1:end-3),1e3*res.(fn{idx4(k)})); %#ok<AGROW>
            end
            axis(out.ax(5),'off');
            hInfo = text(0,1,strInfo,'Parent',out.ax(5),'HorizontalAlignment','left',...
                'VerticalAlignment','top', 'FontName','FixedWidth','Interpreter','none',...
                'FontSize',16);
            % final thoughts
            hLinkAx(1) = linkprop(out.ax([1 3]),'xlim');
            setappdata(out.fig,'hLinkAx',hLinkAx);
            set(out.ax(1:4), 'LineWidth',2,'Fontsize',14);
            set(findall(out.fig,'type','text'),'fontSize',14);
            set(hInfo,'fontSize',16);
            set(out.fig,'Visible','on');
            if nargout > 0
                varargout = {out};
            end
            
            function addFit(str)
                % addFit Perform fit to polynominal based on dX, dY or dR (str) and add parameters
                % according to ISO11146 to output structure
                
                strAdd = @(x) [str x];
                res.(strAdd(''))      = [obj(:).(str)];
                res.(strAdd('_2Fit')) = polyfit(res.zCCD,[obj(:).(str)].^2,2);
                a = res.(strAdd('_2Fit'))(3); b = res.(strAdd('_2Fit'))(2); c = res.(strAdd('_2Fit'))(1);
                res.(strAdd('_z0'))   = -b/2/c;
                res.(strAdd('_d0'))   = 1/(2*sqrt(c)) * sqrt(4*a*c-b^2);
                res.(strAdd('_div'))  = sqrt(c);
                res.(strAdd('_zR'))   = 1/(2*c) * sqrt(4*a*c-b^2);
                res.(strAdd('_M2'))   = pi/(8*mean([obj(:).wavelength])) * sqrt(4*a*c-b^2);
            end
            
            function plotFitD(ax, str)
                strAdd = @(x) [str x];
                plot(ax, 1e3*res.zCCD, 1e3*res.(str), 'Displayname', [str ' - exp'], ...
                    'Linewidth', 2, 'Linestyle','none', 'Marker', 'o', 'Color', mycolor(i,:));
                plot(ax, 1e3*z, 1e3*sqrt(polyval(res.(strAdd('_2Fit')),z)), 'Displayname', [str ' - fit'], ...
                    'Linewidth', 2, 'Linestyle','-', 'Marker', 'none', 'Color', mycolor(i,:));
                i = i + 1;
            end
            
            function plotFitF(ax, str)
                strAdd = @(x) [str x];
                plot(ax, 1e3*res.zCCD, 1e-4*res.(str), 'Displayname', [str ' - exp'], ...
                    'Linewidth', 2, 'Linestyle','none', 'Marker', 'o', 'Color', mycolor(i,:));
                plot(ax, 1e3*z, 1e-4./(polyval(res.(strAdd('_Fit')),z)), 'Displayname', [str ' - fit'], ...
                    'Linewidth', 2, 'Linestyle','-', 'Marker', 'none', 'Color', mycolor(i,:));
                i = i + 1;
            end
            
            function plotWaist(ax)
                
                % plot min, max and mean of:
                %   z0 Location of beam waist
                %   d0 Diameter of beam waist
                %   zR Rayleigh length of waist
                fn  = fieldnames(res);
                idx = cellfun(@(x) numel(x) > 2 && strcmp(x(end-2:end),'_z0'),fn);
                z0  = cellfun(@(x) res.(x),fn(idx));
                idx = cellfun(@(x) numel(x) > 2 && strcmp(x(end-2:end),'_d0'),fn);
                d0  = cellfun(@(x) res.(x),fn(idx));
                idx = cellfun(@(x) numel(x) > 2 && strcmp(x(end-2:end),'_zR'),fn);
                zR  = cellfun(@(x) res.(x),fn(idx));
                % diameter
                plot(ax,1e3*[max(z0-zR/2) min(z0+zR/2)], 1e3*[min(d0) min(d0)],'k--');
                plot(ax,1e3*[mean(z0-zR/2) mean(z0+zR/2)], 1e3*[mean(d0) mean(d0)],'k-');
                plot(ax,1e3*[min(z0-zR/2) max(z0+zR/2)], 1e3*[max(d0) max(d0)],'k--');
                % waist position
                plot(ax,1e3*[min(z0) min(z0)], 1e3*[0 max(d0)],'k--');
                plot(ax,1e3*[mean(z0) mean(z0)], 1e3*[0 max(d0)],'k-');
                plot(ax,1e3*[max(z0) max(z0)], 1e3*[0 max(d0)],'k--');
                % left rayleigh
                plot(ax,1e3*[min(z0-zR/2) min(z0-zR/2)], 1e3*[0 max(d0)],'k--');
                plot(ax,1e3*[mean(z0-zR/2) mean(z0-zR/2)], 1e3*[0 mean(d0)],'k-');
                plot(ax,1e3*[max(z0-zR/2) max(z0-zR/2)], 1e3*[0 min(d0)],'k--');
                % right rayleigh
                plot(ax,1e3*[min(z0+zR/2) min(z0+zR/2)], 1e3*[0 min(d0)],'k--');
                plot(ax,1e3*[mean(z0+zR/2) mean(z0+zR/2)], 1e3*[0 mean(d0)],'k-');
                plot(ax,1e3*[max(z0+zR/2) max(z0+zR/2)], 1e3*[0 max(d0)],'k--');
            end
        end
        
        function              setOrigin(obj,varargin)
            % setOrigin Use a given frame in a given track structure to set CCD position
            % such that the origin is at the center of the track, to set the origin to the center of
            % the beam call the function without any input
            
            if nargin > 1
                % use super class method
                setOrigin@Video(obj,varargin{:});
                return;
            end
            %
            % run one-by-one
            if numel(obj) > 1
                for i = 1:numel(obj)
                    if obj(i).lock
                        warning(sprintf('%s:Input',mfilename),['File ''%s'' is locked to ',...
                            'prevent any data change'],obj(i).filename);
                    else
                        setOrigin(obj(i));
                    end
                end
                return;
            end
            %
            % set origin based on the beam profile
            %
            % determine center of ROI in global coordinate system and set new CCD position such that
            % the origin is at the center of the ROI
            tmp      = [0 0 0];
            tmp([obj.xIdx obj.yIdx]) = [obj.m1X obj.m1Y];
            obj.pos  = obj.pos-tmp;
        end
        
        function         resetUpdate(obj,doNotify)
            % resetUpdate Resets properties that should be recomputed on data change
            %
            % This will force a recalculation of the corresponding properties next time they are
            % used and also reduces the memory consumption of the object
            
            if nargin < 2, doNotify = true; end
            % reset subclass, keep norm2Abs (is supposed to be the same for the complete video) but
            % reset the energy and ISO computation.
            for i = 1:numel(obj)
                obj(i).p_iso11146    = [];
                obj(i).p_energy      = [];
                obj(i).p_ixMax       = [];
                obj(i).p_iyMax       = [];
                obj(i).p_xBeamGrid   = [];
                obj(i).p_yBeamGrid   = [];
                obj(i).p_zBeamGrid   = [];
                obj(i).p_fluenceNorm = [];
                obj(i).p_fluenceMax  = [];
                obj(i).p_fluenceMean = [];
                obj(i).p_fluenceMin  = [];
            end
            % reset superclass
            resetUpdate@Video(obj,doNotify);
        end
    end
    
    methods (Access = protected, Hidden = false)
        function         setProperties(obj,dat,props2Copy)
            %setProperties Sets the properties of the object based on the data in a given structure
            % or another object, but ignores any given fluence data. Use this function to set the
            % properties as returned by readCSV or copy all properties listed in props2Disk from a
            % given object. The third arguments lists the properties or fieldnames to actually copy
            % in case another object is given as source (optional). If 'obj' is scalar but 'dat' is
            % not, the mean value will be computed where possible and set for the property of 'obj'.
            % If 'obj' is not scalar but 'dat' is, all elements of 'obj' are set to the value based
            % on 'dat'.
            
            %
            % check and prepare input
            if ~((isstruct(dat) || isa(dat,class(obj))))
                error(sprintf('%s:Input',mfilename),['Function expects a scalar property structure ',...
                    'or an array structure matching the number of objects, another object instead of ',...
                    'an object will do as well']);
            elseif numel(dat) == 1
                idxObj = 1:numel(obj);
                idxDat = num2cell(ones(size(obj)));
            elseif numel(dat) == numel(obj)
                idxObj = 1:numel(obj);
                idxDat = num2cell(1:numel(obj));
            elseif numel(obj) == 1
                idxObj = 1;
                idxDat = {1:numel(dat)};
            else
                error(sprintf('%s:Input',mfilename),'Mismatch in number of indices, please check!');
            end
            if isa(dat,class(obj))
                % copy properties from another object
                if nargin < 3, props2Copy = obj.props2Disk; end
                for i = 1:numel(idxObj)
                    for n = 1:numel(props2Copy)
                        switch props2Copy{n}
                            case {'name' 'device' 'comment' 'userdata' 'bufferData' 'map' ...
                                    'transform' 'track'}
                                % use data from first object
                                obj(idxObj(i)).(props2Copy{n}) = dat(idxDat{i}(1)).(props2Copy{n});
                            otherwise
                                obj(idxObj(i)).(props2Copy{n}) = mean(cat(1,dat(idxDat{i}).(props2Copy{n})),1);
                        end
                    end
                end
            else
                % set properties from CSV
                for i = 1:numel(idxObj)
                    obj(idxObj(i)).date        = dat(idxDat{i}).date;
                    obj(idxObj(i)).comment     = dat(idxDat{i}).comment;
                    obj(idxObj(i)).pixres      = dat(idxDat{i}).pixres;
                    obj(idxObj(i)).device      = dat(idxDat{i}).device;
                    obj(idxObj(i)).wavelength  = dat(idxDat{i}).wavelength;
                    obj(idxObj(i)).norm2Abs    = mean((10^([dat(idxDat{i}).attenuation]/10))./[dat(idxDat{i}).gain]);
                    myMax = 0;
                    for n = 1:numel(idxDat{i})
                        myMax = max(myMax, double(intmax(class(dat(idxDat{i}(n)).data))));
                    end
                    obj(idxObj(i)).video2Norm  = [mean([dat(idxDat{i}).baseLevel]) myMax];
                end
            end
        end
        
        function value = get_props2Disk(obj)
            %get_props2Disk Return the properties that should be stored to disk, add properties of
            % beam profiles to superclass properties that go to disk
            
            value = get_props2Disk@Video(obj);
            value = unique(cat(2,reshape(value,1,[]),...
                {'video2Norm' 'norm2Abs' 'idxFrames' 'idxChannels' 'wavelength' 'polarization' 'focalLength'}));
        end
        
        function value = get_info(obj)
            % get_info Returns info for a scalar object
            
            %
            % get superclass info
            if numel(obj) > 1
                error(sprintf('%s:Input',mfilename),...
                    'Function only supports scalar inputs');
            end
            value = get_info@Video(obj);
            %
            % add beam specific part
            nIntro   = 24;
            nValue   = 10;
            props    = {'idxFrames' 'idxChannels' 'energy' 'wavelength' 'polarization' 'focalLength' 'sep' ...
                'm1X' 'dX' 'AFlatMax' 'm2X' 'xBeamMax' 'sep' ...
                'rot' 'ell' 'ecc' 'sep' ...
                'fluenceMax' 'fluenceMin' 'fluenceMean'
                };
            value{end+1} = sprintf('%*s',round(1.5*nIntro),'Properties of beam profile');
            for i = 1:numel(props)
                switch props{i}
                    case 'sep'
                        value{end+1} = ''; %#ok<AGROW>
                    case 'idxFrames'
                        if numel(obj.idxFrames) <= 10
                            str = sprintf('[%s]',num2str(obj.idxFrames));
                        else
                            str = sprintf('%d frames, range: %d - %d',...
                                numel(obj.idxFrames),min(obj.idxFrames),max(obj.idxFrames));
                        end
                        value{end+1} = sprintf('%*s: %s',nIntro,'selected frame(s)',str); %#ok<AGROW>
                    case 'idxChannels'
                        if numel(obj.idxChannels) <= 10
                            str = sprintf('[%s]',num2str(obj.idxChannels));
                        else
                            str = sprintf('%d channels, range: %d - %d',...
                                numel(obj.idxChannels),min(obj.idxChannels),max(obj.idxChannels));
                        end
                        value{end+1} = sprintf('%*s: %s',nIntro,'selected channel(s)',str); %#ok<AGROW>
                    case 'energy'
                        value{end+1} = sprintf('%*s: %*.2e mJ',nIntro,props{i},nValue,1e3*obj.(props{i})); %#ok<AGROW>
                    case 'wavelength'
                        value{end+1} = sprintf('%*s: %*.2f nm',nIntro,props{i},nValue,1e9*obj.(props{i})); %#ok<AGROW>
                    case 'focalLength'
                        value{end+1} = sprintf('%*s: %*.2f mm',nIntro,props{i},nValue,1e3*obj.(props{i})); %#ok<AGROW>
                    case 'm1X'
                        str          = ['[' num2str([obj.m1X obj.m1Y]*1e3) '] mm'];
                        value{end+1} = sprintf('%*s: %s',nIntro,'position [m1X m1Y  ]',str); %#ok<AGROW>
                    case 'm2X'
                        str          = ['[' num2str([obj.m2X obj.m2Y obj.m2XY]*1e6) '] mm^2'];
                        value{end+1} = sprintf('%*s: %s',nIntro,'moments [m2X m2Y m2XY]',str); %#ok<AGROW>
                    case 'dX'
                        str          = ['[' num2str([obj.dX obj.dY]*1e3) '] mm'];
                        value{end+1} = sprintf('%*s: %s',nIntro,'size  [dX  dY  ]',str); %#ok<AGROW>
                        str          = ['[' num2str([obj.dR obj.dEff]*1e3) '] mm'];
                        value{end+1} = sprintf('%*s: %s',nIntro,'[dR  dEff]',str); %#ok<AGROW>
                        str          = ['[' num2str([obj.dFlatMax obj.dFlatMean]*1e3) '] mm'];
                        value{end+1} = sprintf('%*s: %s',nIntro,'[dFlatMax dFlatMean]',str); %#ok<AGROW>
                    case 'AFlatMax'
                        str          = ['[' num2str([obj.AFlatMax obj.AFlatMean]*1e6) '] mm^2'];
                        value{end+1} = sprintf('%*s: %s',nIntro,'[AFlatMax AFlatMean]',str); %#ok<AGROW>
                    case 'xBeamMax'
                        str          = ['[' num2str([obj.xBeamMax obj.yBeamMax]*1e3) '] mm^2'];
                        value{end+1} = sprintf('%*s: %s',nIntro,'max [xBeamMax yBeamMax]',str); %#ok<AGROW>
                    case {'rot' 'polarization'}
                        value{end+1} = sprintf(['%*s: %*.2f ' char(176)],nIntro,props{i},nValue,180/pi*obj.(props{i})); %#ok<AGROW>
                    case {'ell' 'ecc'}
                        value{end+1} = sprintf('%*s: %*.2g ',nIntro,props{i},nValue, obj.(props{i})); %#ok<AGROW>
                    case {'fluenceMax' 'fluenceMin' 'fluenceMean'}
                        value{end+1} = sprintf('%*s: %*.2e J/cm^2',nIntro,props{i},nValue,1e-4*obj.(props{i})); %#ok<AGROW>
                    otherwise
                        if isnumeric(obj.(props{i})) && isscalar(obj.(props{i}))
                            value{end+1} = sprintf('%*s: %*.2e',nIntro,props{i},nValue,obj.(props{i})); %#ok<AGROW>
                        elseif ischar(obj.(props{i}))
                            value{end+1} = sprintf('%*s: %s',nIntro,props{i},obj.(props{i})); %#ok<AGROW>
                        elseif iscellstr(obj.(props{i}))
                            value{end+1} = sprintf('%*s: %s',nIntro,props{i},obj.(props{i}){1}); %#ok<AGROW>
                            for j = 2:numel(obj.(props{i}))
                                value{end+1} = sprintf('%*s  %s',nIntro,' ',obj.(props{i}){j}); %#ok<AGROW>
                            end
                        end
                end
            end
            value{end+1} = value{1};
        end
        
        function value = get_memory(obj)
            % get_memory Returns memory for a scalar object
            
            if numel(obj) > 1
                error(sprintf('%s:Input',mfilename),...
                    'Function only supports scalar inputs');
            end
            % memory in supclass
            value = get_memory@Video(obj);
            % add data from subclass
            value = value + ...
                sum(double([~isempty(obj.p_fluenceNorm) ~isempty(obj.p_xBeamGrid) ...                       % fluence and grid
                ~isempty(obj.p_yBeamGrid) ~isempty(obj.p_zBeamGrid)])) * obj.nX * obj.nY * 8 / 1024^2 + ... % more grids
                (1+double(~isempty(obj.p_backupData))) * 100 * 8 / 1024^2;                                  % estimate all other properties
        end
    end
    
    %% Static class related methods
    methods (Static = true, Access = public, Hidden = false)
        function obj       = loadobj(S)
            % loadobj Loads object and recalls data, necessary to load beam profiles that were not a
            % subclass of Video
            
            %
            % call supclass in case it is a regular object
            if isstruct(S) && ~(isfield(S,'fluenceNorm') || isfield(S,'p_fluenceNorm'))
                obj = loadobj@Video(S);
                return
            end
            %
            % the load process should be accessing an older version, that needs to be converted now
            % to the new class, this is done by creating a BeamProfile object without memory mapping
            fn = fieldnames(S);
            if all(ismember({'p_norm2Abs' 'p_device' 'p_comment' 'p_fluenceNorm' 'p_pos' 'p_energy',...
                    'p_time' 'p_pixres' 'p_lambda' 'p_pol' 'p_focalLength' 'p_userdata'},fn))
                nObj = numel(S);
                obj  = BeamProfile(nObj);
                for n = 1:nObj
                    obj(n).cdata.cdata  = S(n).p_fluenceNorm;
                    obj(n).video2Norm   = [0 1];
                    obj(n).idxFrames    = 1;
                    obj(n).idxChannels  = 1;
                    obj(n).name         = 'Loaded beam profile written by older class definition';
                    obj(n).norm         = [0 0 -1];
                    obj(n).device       = S(n).p_device;
                    obj(n).comment      = S(n).p_comment;
                    obj(n).userdata     = S(n).p_userdata;
                    obj(n).date         = S(n).p_time;
                    if ~isempty(S(n).p_norm2Abs)
                        obj(n).norm2Abs = S(n).p_norm2Abs;
                    end
                    if ~isempty(S(n).p_pos) && all(~isnan(S(n).p_pos))
                        obj(n).pos = S(n).p_pos;
                    end
                    if ~isempty(S(n).p_pixres) && all(~isnan(S(n).p_pixres))
                        obj(n).pixres = S(n).p_pixres;
                    end
                    obj(n).energy       = S(n).p_energy;
                    obj(n).wavelength   = S(n).p_lambda;
                    obj(n).polarization = S(n).p_pol;
                    obj(n).focalLength  = S(n).p_focalLength;
                end
            elseif all(ismember({'device' 'comment' 'fluenceNorm' 'pos' 'energy',...
                    'time' 'pixres' 'lambda' 'focalLength' 'userdata'},fn))
                nObj = numel(S);
                obj  = BeamProfile(nObj);
                for n = 1:nObj
                    obj(n).cdata.cdata  = S(n).fluenceNorm;
                    obj(n).video2Norm   = [0 1];
                    obj(n).idxFrames    = 1;
                    obj(n).idxChannels  = 1;
                    obj(n).name         = 'Loaded beam profile written by older class definition';
                    obj(n).norm         = [0 0 -1];
                    obj(n).device       = S(n).device;
                    obj(n).comment      = S(n).comment;
                    obj(n).userdata     = S(n).userdata;
                    obj(n).date         = S(n).time;
                    if isfield(S,'p_norm2Abs') && ~isempty(S(n).p_norm2Abs)
                        obj(n).norm2Abs = S(n).p_norm2Abs;
                    end
                    if isfield(S,'pol') && ~isempty(S(n).pol)
                        obj(n).polarization = S(n).pol;
                    end
                    if ~isempty(S(n).pos) && all(~isnan(S(n).pos))
                        obj(n).pos = S(n).pos;
                    end
                    if ~isempty(S(n).pixres) && all(~isnan(S(n).pixres))
                        obj(n).pixres = S(n).pixres;
                    end
                    obj(n).energy       = S(n).energy;
                    obj(n).wavelength   = S(n).lambda;
                    obj(n).focalLength  = S(n).focalLength;
                end
            else
                error(sprintf('%s:Input',mfilename),['Unknown structure during load process of a ',...
                    'BeamProfile object, maybe even older data that was not considered so far?']);
            end
        end
        
        function varargout = placeDropPlot(out,varargin)
            % placeDropPlot Plot the result of the placeDrop function
            
            nargoutchk(0,2);
            %
            % check and parse input
            if ~(isstruct(out) && isscalar(out))
                error(sprintf('%s:Input',mfilename),['Expected a scalar structure with the results ',...
                    'of placeDrop as input, please check!']);
            end
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % property to plot finally
            opt.addParameter('plotme', 'convCFAbs', ...
                @(x) ischar(x));
            % true/false whether to normalize the axi-symmetric result of this function (not the
            % plotted profile) to the cos * gaussian fit, that way the results of an actual and
            % prefect beam should be easier to compare
            opt.addParameter('plotNorm', true, ...
                @(x) islogical(x) && isscalar(x));
            % number of points to plot some objects
            opt.addParameter('nPlot', 100, ...
                @(x) isnumeric(x) && isscalar(x));
            % true/false whether to plot beam profile with the drop
            opt.addParameter('plotProfile', true, ...
                @(x) islogical(x) && isscalar(x));
            % true/false whether to reuse existing figure with correct tag(s)
            opt.addParameter('plotReuse', true, ...
                @(x) islogical(x) && isscalar(x));
            opt.parse(varargin{:});
            opt = opt.Results;
            %
            % plot results
            if ~(ismember(opt.plotme,fieldnames(out)) && ismatrix(out.(opt.plotme)))
                error(sprintf('%s:Input',mfilename),'Unknown property ''%s'' to plot, please check!',opt.plotme);
            end
            if numel(out.parCOSGauss) > 2
                myModel = @BeamProfile.myCOSGaussianNonParallel;
            else
                myModel = @BeamProfile.myCOSGaussian;
            end
            %
            % normalize normalized fluence with cos*gaussian at theta = 0
            if opt.plotNorm
                fNorm = 1/myModel(out.parCOSGauss,0);
            else
                fNorm = 1;
            end
            % plot axi-symmtric representation
            out.fig = [];
            if opt.plotReuse
                fig = findall(groot,'type','figure', 'Tag','BeamProfile_placeDrop');
                if ~isempty(fig) && isvalid(fig(1))
                    out.fig = fig(1);
                    delete(out.fig.Children);
                end
            end
            if isempty(out.fig)
                out.fig = figure('Name',sprintf('Drop (R0=%.4f mm) in beam profile from %s', ...
                    1e3*out.R0, datestr(out.beamprofile.date)), 'Visible','off',...
                    'Tag','BeamProfile_placeDrop');
            end
            out.ax  = subplot(2,3,[1 4],'Nextplot','Add','Parent',out.fig);
            phiPlot = linspace(0,out.phiMax,100);
            plot(out.ax(1), 180/pi*out.linPhi, out.linCF*fNorm,'Displayname','exp - all','Color',0.8*[1 1 1]);
            plot(out.ax(1), 180/pi*smooth(out.linPhi,10), smooth(out.linCF,10)*fNorm,'Displayname','exp - smooth 5','Color',0.6*[1 1 1]);
            plot(out.ax(1), 180/pi*smooth(out.linPhi,50), smooth(out.linCF,50)*fNorm,'Displayname','exp - smooth 10','Color',0.3*[1 1 1]);
            plot(out.ax(1), 180/pi*out.fitPhi, out.fitCF*fNorm,'xr','Displayname','exp - mean');
            plot(out.ax(1), 180/pi*phiPlot, BeamProfile.myGaussian(out.parGauss,phiPlot)*fNorm,...
                '--','Displayname',sprintf('gaussian, sig = %.2f pi',out.parGauss(2)/pi),'Linewidth',2);
            plot(out.ax(1), 180/pi*phiPlot, myModel(out.parCOSGauss,phiPlot)*fNorm,...
                '-','Displayname',sprintf('cos * gaussian, sig = %.2f pi',out.parCOSGauss(2)/pi),'Linewidth',2);
            xlabel(out.ax(1),'rot (deg)');
            ylabel(out.ax(1),'norm. fluence');
            xlim(out.ax(1),[0 out.phiMax*180/pi]);
            legend(out.ax(1),'show');
            % plot gridded data, e.g. the transmitted fluence
            idxI = max(1,find(sum(out.mask,2),1,'first')):min(size(out.mask,1),find(sum(out.mask,2),1,'last'));
            idxJ = max(1,find(sum(out.mask,1),1,'first')):min(size(out.mask,2),find(sum(out.mask,1),1,'last'));
            out.ax(2) = subplot(2,3,[ 2 3 5 6],'Nextplot','Add','Parent',out.fig);
            imagesc(out.beamprofile.xDir * [-1 1]*out.rMax*1e3,out.beamprofile.yDir * [-1 1]*out.rMax*1e3,out.(opt.plotme)(idxI,idxJ),...
                'Parent',out.ax(2));
            % add drop and aperture
            xE = @(t,x,y,r) x + r * cos(t);
            yE = @(t,x,y,r) y + r * sin(t);
            t  = linspace(0,2*pi,opt.nPlot);
            plot(out.ax(2),xE(t,0,0,out.R0)*1e3,yE(t,0,0,out.R0)*1e3,...
                'Color','black','Linestyle','--','DisplayName','Drop');
            plot(out.ax(2),xE(t,0,0,out.rMax)*1e3,yE(t,0,0,out.rMax)*1e3,...
                'Color','black','Linestyle','-','DisplayName','Drop aperture');
            colorbar(out.ax(2),'northoutside');
            axis(out.ax(2), 'equal');
            strX = ['drop x, global ' out.beamprofile.xStr ' (mm)'];
            strY = ['drop y, global ' out.beamprofile.yStr ' (mm)'];
            xlabel(out.ax(2),strX);
            ylabel(out.ax(2),strY);
            title(out.ax(2),opt.plotme);
            if out.beamprofile.xDir > 0, set(out.ax(2),'XDir','normal');
            else,                        set(out.ax(2),'XDir','reverse');
            end
            if out.beamprofile.yDir > 0, set(out.ax(2),'YDir','reverse');
            else,                        set(out.ax(2),'YDir','normal');
            end
            set(out.fig,'Visible','on','HandleVisibility','callback');
            out.fig = reshape(out.fig,1,[]);
            out.ax  = reshape(out.ax,1,[]);
            %
            % plot drop position on beam profile
            if opt.plotProfile
                hGraph = out.beamprofile.plotProfile('plotReuse',opt.plotReuse);
                xE     = @(t,x,y,r) x + r * cos(t);
                yE     = @(t,x,y,r) y + r * sin(t);
                t      = linspace(0,2*pi,opt.nPlot);
                plot(hGraph.ax(1), xE(t,out.beamprofile.m1X+out.xy0(1),out.beamprofile.m1Y+out.xy0(2),out.R0)*1e3, ...
                    yE(t,out.beamprofile.m1X+out.xy0(1),out.beamprofile.m1Y+out.xy0(2),out.R0)*1e3, ...
                    'Color','black','Linestyle','--','DisplayName','Drop');
                plot(hGraph.ax(1), xE(t,out.beamprofile.m1X+out.xy0(1),out.beamprofile.m1Y+out.xy0(2),out.rMax)*1e3, ...
                    yE(t,out.beamprofile.m1X+out.xy0(1),out.beamprofile.m1Y+out.xy0(2),out.rMax)*1e3, ...
                    'Color','black','Linestyle','-','DisplayName','Drop aperture');
                if isfield(out,'fig')
                    out.fig = cat(2,out.fig,reshape(hGraph.fig,1,[]));
                    out.ax   =cat(2,out.ax,reshape(hGraph.ax,1,[]));
                else
                    out.fig = reshape(hGraph.fig,1,[]);
                    out.ax  = reshape(hGraph.ax,1,[]);
                end
            end
            if nargout > 0
                varargout = {out.fig out.ax};
                varargout = varargout(1:nargout);
            end
        end
        
        function ydata     = myGaussian(par,xdata)
            %myGauss Gaussian with two parameters: amplitude and std
            ydata = par(1) * exp(-xdata.^2/(2*par(2)^2));
        end
        
        function ydata     = myCOSGaussian(par,xdata)
            %myCOSGaussian Cos convoluted with a gaussian with two parameters: amplitude and std
            ydata = par(1) * cos(xdata) .* exp(-xdata.^2/(2*par(2)^2));
        end
        
        function ydata     = myCOSGaussianNonParallel(par,xdata)
            %myCOSGaussianNonParallel Cos convoluted with a gaussian with three parameters: amplitude, std
            % and z0/R0 (the drop position normalized with its initial radius)
            ydata = par(1) * ...
                cos(xdata + atan(sin(xdata)./(sign(par(3))*(abs(par(3))+cos(xdata))))) ... % cosine effect
                .* exp(-xdata.^2/(2*par(2)^2));       % spatial distribution in beam profile
        end
        
        function varargout = convertCSV2VID(filename,varargin)
            % convertCSV2VID Convert given CSV files to a single video taking all beam related
            % properties (e.g. pixel resolution, wavelength, etc.) from the first file and the
            % actual CCD data from all files.
            %
            % For available input options have a look at the actual workhorse function,
            % Videomap.convert2VID. This function justs sets the import option, such that
            % convert2VID can process CSV files.
            %
            % The output returned by this function is the filename of the new file, please note: the
            % beam profile information is stored in a video file (e.g. TIF) and a MAT file that
            % holds all additional data except the CCD image data.
            
            if nargin < 1, varargout = {}; return; end
            nargoutchk(0,1);
            %
            % capture only the import option and let everything pass to conert2Vid in Videomap class
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = true;
            % a custom read function to read given files, should return a single image when a single
            % filename is given as input,
            opt.addParameter('import', [], ...
                @(x) isempty(x) || isa(x,'function_handle'));
            opt.parse(varargin{:});
            tmp = opt.Unmatched;
            fn  = fieldnames(opt.Results);
            for n = 1:numel(fn)
                tmp.(fn{n}) = opt.Results.(fn{n});
            end
            opt = tmp;
            if ~isempty(opt.import)
                warning(sprintf('%s:Conversion',mfilename),['Setting an import function is not expected, ',...
                    ' since the convertCSV2VID function takes care of it normally, please be careful!']);
            else
                opt.import = @myread;
            end
            %
            % convert to a video file, load it and set beam specific options based on the first CSV
            % file that was read in convert2VID
            [out, fnCSV] = BeamProfile.convert2VID(filename, opt);
            bp           = BeamProfile(out);
            prop         = BeamProfile.readCSV(fnCSV(1).name);
            setProperties(bp,prop);
            bp.store;
            fprintf(['  Beam profile properties set based on file''%s''\n',...
                '  Beam profile data stored in ''%s.mat''\n'],fnCSV(1).name,bp.filename);
            if nargout > 0, varargout= {out}; end
            
            function img = myread(fn)
                dat = BeamProfile.readCSV(fn);
                img = dat.data;
            end
        end
    end
    
    methods (Access = protected)
        function cpObj = copyElement(obj)
            % copyElement Override copyElement method from matlab.mixin.Copyable class
            
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % reset new object
            cpObj.resetUpdate;
        end
    end
    
    %% Static class related methods
    methods (Static = true, Access = public, Hidden = false, Sealed = true)
        function dat = readRAW(filename)
            %readRAW Reads raw image from Thorlabs stored in 12 bit mode
            %   This is a reverse engineered read function without much error checking, please be
            %   careful when using it, it also returns just the image data and no additional info
            
            narginchk(1,1);
            assert(ischar(filename) && exist(filename,'file') ==2,...
                'Could not find data file with given input, please check!');
            fn   = dir(filename);
            fn   = fn(1);
            info = imfinfo(filename);
            w    = abs(info.Width);
            h    = abs(info.Height);
            iS   = fn.bytes/2-w*h+1;
            % check if file has reached a minimum size and assume the image data is in the last bytes
            if fn.bytes < 1024 || iS < 1
               counter = 0; maxcount = 20; doAgain = true; delay = 0.2; pause(delay);
               while doAgain
                   fn      = dir(filename);
                   fn      = fn(1);
                   if fn.bytes < 1024
                       doAgain = true;
                   else
                       info    = imfinfo(filename);
                       w       = abs(info.Width);
                       h       = abs(info.Height);
                       iS      = fn.bytes/2-w*h+1;
                       doAgain = iS < 1;
                   end
                   doAgain = doAgain && counter < maxcount;
                   pause(delay);
                   counter = counter + 1;
               end
               if iS < 1
                   error(sprintf('%s:Input',mfilename),'Data file ''%s'' did not reached minimum size',filename);
               end
            end
            % read data and return last bytes
            fid  = fopen(filename);
            dat  = fread(fid,'uint16=>uint16');
            fclose(fid);
            dat = (reshape(dat(iS:(iS+h*w-1)),[w h])');
        end
        
        function dat = readCSV(filename)
            %readCSV Reads CCD data from CSV file written by Thorlabs Beam software
            %
            % Assumptions for import:
            %   * Use importdata to read file with header and data
            %   * Normalize fluence from CCD with maximum value and store total energy in the beam
            %   * Pixel resolution according to manual for 'BC106-VIS' or 'BC106N-VIS/M'
            %
            % Output structure:
            %   attenuation Attenuation in db (double)
            %   comment     Comments in file (cellstr)
            %   baseLevel   Base level (double)
            %   data        CCD data as uint8 or unit16 aligned to MSB (uint8 or uint16)
            %   date        Time of measurement (datetime)
            %   device  	Name of beam profiler and serial number (string)
            %   exposure    Exposure of CCD in s (double)
            %   gain        Gain of CCD (double)
            %   nBits       Number of bits per pixel during recording (double)
            %   pixres      Pixel resolution of CCD chip in m/pix in horizontal and vertical direction (2 x double)
            %   roi         Region of interest as a MATAB rect (double)
            %   sizeCCD     CCD size in pixel in horizontal and vertival direction (2 x double)
            %   wavelength  Wavelength of beam in nm (double)
            
            %
            % check input
            narginchk(1,1);
            [~,~,ext] = fileparts(filename);
            if exist(filename,'file') ~= 2
                error(sprintf('%s:Input',mfilename),'Could not find data file ''%s''',filename);
            elseif ~strcmp(ext,'.csv')
                error(sprintf('%s:Input',mfilename),'Data file ''%s'' is not a CSV file',filename);
            end
            %
            % read file
            tmp         = importdata(filename);
            out.Comment = {};
            for i = 1:numel(tmp.textdata)
                str = strtrim(tmp.textdata{i});
                if isempty(str)
                    % do nothing
                elseif strcmp(str(1),'%') || strcmp(str(1),'#')
                    % add comment
                    out.Comment{end+1} = str;
                else
                    [token, remain] = strtok(str,':');
                    token  = strtrim(token);
                    if ~isempty(token)
                        remain = strtrim(remain(2:end));
                    end
                    out.(matlab.lang.makeValidName(token)) = remain;
                end
            end
            %
            % check content, read mandotory values and assign defaults for remaining properties
            if ~all(ismember({'ThorlabsBeam' 'Comment' 'Version' 'Date' 'Time' 'Device' 'S_N'},fieldnames(out)))
                error(sprintf('%s:Input',mfilename),'File ''%s'' does not seem to contain all required information',filename);
            end
            dat.comment     = cat(1,{sprintf('File ''%s'' written with Thorlabs Beam (version: %s)',filename,out.Version)}, out.Comment{:});
            dat.date        = datetime(sprintf('%s %s',out.Date,out.Time),'InputFormat','M/d/uuuu h:m:s a');
            dat.device      = sprintf('%s - S/N:%s',out.Device,out.S_N);
            dat.pixres      = NaN(1,2);
            dat.wavelength  = NaN;
            dat.attenuation = NaN;
            dat.exposure    = NaN;
            dat.gain        = NaN;
            dat.nBits       = NaN;
            dat.baseLevel   = NaN;
            dat.roi         = NaN(1,4);
            dat.sizeCCD     = NaN(1,2);
            fn              = fieldnames(out);
            %
            % parse text
            for i = 1:numel(fn)
                switch fn{i}
                    case 'PixelSize__m_'
                        res = regexp(out.(fn{i}),'Horizontal: (?<pixresH>\d+(\.\d{1,2})?), Vertical: (?<pixresV>\d+(\.\d{1,2})?)','names');
                        if ~isempty(res)
                            dat.pixres = [str2double(res.pixresH) str2double(res.pixresV)]*1e-6;
                        end
                    case 'Resolution'
                        res = regexp(out.(fn{i}),'(?<nbits>\d+) bit','names');
                        if ~isempty(res)
                            dat.nBits = str2double(res.nbits);
                        end
                    case 'Wavelength_nm_'
                        dat.wavelength = str2double(out.(fn{i})) * 1e-9;
                    case 'Attenuation_dB_'
                        dat.attenuation = str2double(out.(fn{i}));
                    case 'ExposureTime_ms_'
                        dat.exposure = str2double(out.(fn{i}))*1e-3;
                    case 'Gain'
                        dat.gain = str2double(out.(fn{i}));
                    case 'BaseLevel_digits_'
                        dat.baseLevel = str2double(out.(fn{i}));
                    case 'RegionOfInterest_pix_'
                        res = regexp(out.(fn{i}),'Horizontal: (?<r1>\d+), Vertical: (?<r2>\d+), Width (?<w>\d+), Height (?<h>\d+)','names');
                        if ~isempty(res)
                            dat.roi = [str2double(res.r1)+1 str2double(res.r2)+1 str2double(res.w) str2double(res.h)];
                        end
                    case 'SensorResolution_pix_'
                        res = regexp(out.(fn{i}),'Width: (?<w>\d+), Height: (?<h>\d+)','names');
                        if ~isempty(res)
                            dat.sizeCCD = [str2double(res.w) str2double(res.h)];
                        end
                end
            end
            % fallback if information is missing in files
            if any(isnan(dat.pixres))
                switch out.Device
                    case {'BC106-VIS' 'BC106N-VIS/M'}
                        warning(sprintf('%s:Input',mfilename),['Data file ''%s'' seems to lack some information, ',...
                            'but the device ''%s'' was recognized, please check nevertheless!'],filename,out.Device);
                        dat.pixres = [1 1] * 6.45e-6;
                        if any(isnan(dat.sizeCCD))
                            dat.sizeCCD = [1360 1024];
                        end
                    otherwise
                        error(sprintf('%s:Input',mfilename),['Data file ''%s'' seems to lack some information ',...
                            'and the device ''%s'' was NOT recognized, please check!'],filename,out.Device);
                end
            end
            %
            % check data
            if all(~isnan(dat.roi))
                siz = size(tmp.data);
                if ~isequal(siz(1:2),dat.roi([4 3]))
                    warning(sprintf('%s:Input',mfilename),['Size of data ([%s]) in file ''%s'' does not match, ',...
                        'the expected size ([%s]), please check!'], num2str(siz(1:2)), filename, dat.roi([4 3]));
                end
            end
            %
            % convert data
            dataOK = true;
            switch dat.nBits
                case 8
                    dataOK   = max(tmp.data) <= intmax('uint8');
                    dat.data = uint8(tmp.data);
                case 12
                    dataOK   = max(tmp.data) < 2^12;
                    dat.data = bitshift(uint16(tmp.data),4);
                case 16
                    dataOK   = max(tmp.data) <= intmax('uint16');
                    dat.data = uint16(tmp.data);
                otherwise
                    warning(sprintf('%s:Input',mfilename),['Bit resolution for file ''%s'' is unknown, ',...
                        'converting to uint16 without any bitshift, trying to estimate bit resolution, please check'],filename);
                    mymax = max(tmp.data(:));
                    if mymax > intmax('uint16')
                        dat.data = double(tmp.data);
                    elseif mymax > intmax('uint8')
                        dat.data = uint16(tmp.data);
                    else
                        dat.data = uint8(tmp.data);
                    end
            end
            if ~dataOK
                warning(sprintf('%s:Input',mfilename),['Data in file ''%s'' seems to be out of bound ',...
                    'given the bit resolution, please check'],filename);
            end
            dat = orderfields(dat);
        end
        
        function dat       = readMAT(filename)
            %readMAT Reads a single beam profile from a MAT file
            %
            % Storing beam profiles to MAT files is not used directly anymore by this class, since
            % it is now a subclass of the Video class, that takes care of storing the data. However
            % this function can be used to load old files and recover data.
            %
            % Assumptions for import:
            %   * MAT file contains one field with a structure or BeamProfile object OR
            %   * MAT file contains multiple variables but the beam profile is called BeamProfile
            %
            % Output structure:
            %   depends on the MAT file
            %
            
            %
            % check input
            narginchk(1,1);
            [~,~,ext] = fileparts(filename);
            if exist(filename,'file') ~= 2
                error(sprintf('%s:Input',mfilename),'Could not find data file ''%s''',filename);
            elseif ~strcmp(ext,'.mat')
                error(sprintf('%s:Input',mfilename),'Data file ''%s'' is not a MAT file',filename);
            end
            %
            % read file and check if it is either a BeamProfile object or a structure with fields
            % than are understood by BeamProfile class
            tmp          = load(filename);
            fn           = fieldnames(tmp);
            [idxOK, idx] = ismember('BeamProfile',fn);
            if numel(fn) == 1 && (isa(tmp.(fn{1}),'BeamProfile') || ...
                    (isstruct(tmp.(fn{1})) && all(ismember({'device','fluenceNorm','pixres','energy'},fieldnames(tmp.(fn{1}))))))
                dat = tmp.(fn{1});
            elseif numel(fn) > 1 && idxOK && (isa(tmp.(fn{idx}),'BeamProfile') || ...
                    (isstruct(tmp.(fn{idx})) && all(ismember({'device','fluenceNorm','pixres','energy'},fieldnames(tmp.(fn{idx}))))))
                dat = tmp.(fn{idx});
            elseif numel(fn) == 1 && ...
                    (isstruct(tmp.(fn{1})) && all(ismember({'device','data','pixres'},fieldnames(tmp.(fn{1})))))
                warning(sprintf('%s:Input',mfilename),...
                    'File ''%s'' does not include the total energy and may include multiple beam profiles (only one is read)',filename);
                dat.device      = tmp.(fn{1}).device;
                dat.fluenceNorm = tmp.(fn{1}).data(:,:,1);
                dat.fluenceNorm = dat.fluenceNorm./max(dat.fluenceNorm(:));
                dat.pixres      = tmp.(fn{1}).pixres;
            else
                error(sprintf('%s:Input',mfilename),'File ''%s'' does not seem to contain all required information',filename);
            end
            if isstruct(dat)
                dat = orderfields(dat);
            end
        end
        
        function dat       = HermiteGaussianBeam(varargin)
            % HermiteGaussianBeam Creates 3D fluence based on Hermite-Gaussian modes, use
            % implay(abs(dat)) to follow the beam in an animation, it is possible to add several
            % modes by giving vectors for m and n,
            %
            % Needs more checking, x and y axes correct?
            
            % parse and check input
            opt               = inputParser;
            opt.StructExpand  = true;
            opt.KeepUnmatched = false;
            % mode in x
            opt.addParameter('m', 0, ...
                @(x) isnumeric(x) && isvector(x));
            % mode in y
            opt.addParameter('n', 0, ...
                @(x) isnumeric(x) && isvector(x));
            % number of pixels in x
            opt.addParameter('nX', 1024, ...
                @(x) isnumeric(x) && isscalar(x));
            % number of pixels in y
            opt.addParameter('nY', 1360, ...
                @(x) isnumeric(x) && isscalar(x));
            % pixel resolution in m/pix
            opt.addParameter('pixres', 6.45e-6, ...
                @(x) isnumeric(x) && isscalar(x));
            % point on z axis to compute fluence in m
            opt.addParameter('z', 0, ...
                @(x) isnumeric(x) && isvector(x));
            % beam waist in x in m
            opt.addParameter('w0X', 20e-6, ...
                @(x) isnumeric(x) && isscalar(x));
            % beam waist in y in m
            opt.addParameter('w0Y', 20e-6, ...
                @(x) isnumeric(x) && isscalar(x));
            % wavelength in m
            opt.addParameter('wavelength', 532e-9, ...
                @(x) isnumeric(x) && isscalar(x));
            % return normalized absolute part
            opt.addParameter('doAbsNorm', false, ...
                @(x) islogical(x) && isscalar(x));
            opt.parse(varargin{:});
            opt = opt.Results;
            if numel(opt.n) ~= numel(opt.m)
                error(sprintf('%s:Input',mfilename),'m and n must be the same length');
            end
            % create fluence
            x            = ((1:opt.nX)-(1+opt.nX)/2)*opt.pixres;
            y            = ((1:opt.nY)-(1+opt.nY)/2)*opt.pixres;
            [xg, yg, zg] = ndgrid(x,y,opt.z);
            dat = un(xg,zg,opt.m(1),opt.w0X,opt.w0X^2*pi/opt.wavelength).*un(yg,zg,opt.n(1),opt.w0Y,opt.w0Y^2*pi/opt.wavelength);
            for i = 2:numel(opt.n)
                dat = dat + un(xg,zg,opt.m(i),opt.w0X,opt.w0X^2*pi/opt.wavelength).*un(yg,zg,opt.n(i),opt.w0Y,opt.w0Y^2*pi/opt.wavelength);
            end
            if opt.doAbsNorm
                dat = abs(dat);
                dat = dat/max(dat(:));
            end
            dat = permute(dat,[2 1 3]);
            
            function u = un(x,z,n,w0,zR)
                q  = @(z) z + 1i*zR;
                qc = @(z) z - 1i*zR;
                w  = @(z) w0*sqrt(1+(z/zR).^2);
                u  = (2/pi)^0.25 * (2^n*factorial(n)*w0)^-0.5 * (q(0)./q(z)).^0.5 .* ...
                    (q(0)/qc(0) * qc(z)./q(z)).^(n/2) .* BeamProfile.hermite(n,sqrt(2)*x./w(z)) .* ...
                    exp(-1i*zR/w0^2*x.^2./q(z));
            end
        end
        
        function y         = hermite(n,x)
            % hermite Returns the Nth degree Hermite polynomial (physicists' form) of X, where N is
            % a non-negative integer and X a real scalar.
            
            if ~(isnumeric(n) && isscalar(n) && n>=0)
                error(sprintf('%s:Input',mfilename),'n should be a non-negative scalar integer.');
            elseif ~isnumeric(x)
                error(sprintf('%s:Input',mfilename),'x should be numeric.');
            end;
            n = round(n);
            if n == 0
                y = ones(size(x));
            elseif n==1
                y = 2*x;
            elseif n > 1
                y = 2*x.*BeamProfile.hermite(n-1,x) - 2*(n-1)*BeamProfile.hermite(n-2,x);
            end;
        end
    end
end
