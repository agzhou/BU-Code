% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% VSX (Verasonics Script eXecute) for processing .mat file (filename), for
% use exclusively on Verasonics "Vantage" system products.  This version of
% VSX is not compatible with the older "V1" product line.
%
% Usage: run VSX and input 'filename'
%   where filename is the name of a .mat file containing all the structures
%   needed.
%
%   To allow automated scripts to launch VSX, such as HardwareTest scripts,
%   'filename' may be set in the workspace before calling VSX.  To prevent
%   subsequently run scripts from inadvertently using that same 'filename',
%   VSX will use the detected 'filename' once and then delete the variable.
%
%
% The keyword "TEST" has been inserted in locations where temporary edits
% can be made to facilitate debug or testing.
%
% Program steps;
%   1. Read in structures needed to set up for running sequence.
%      - Use Resource structure to set up RcvData, ImgData, and
%        IQData (if needed.)
%   2. Validate structures to see if consistent.
%   3. Add VDAS parameters if missing and running on VDAS hardware.
%   3. Open GUI window.
%   4. Set up figure window for output display.
%   5. Open VDAS hardware, if available. If hardware not available, use
%      simulate mode. To skip check for hardware, set Resource.Parameters.simulateMode=1.
%   6. Enter processing loop:
%      While exit==0
%         Call runAcq to process events in sequence.
%         Every jump command back to 1st event will return control to Matlab
%            to check for GUI actions.
%            - for Freeze action, stop sequencer and wait in Matlab for
%              Freeze button to change state.
%   7. On exit, copy back processing buffers from C to Matlab.
%   8. Close hardware if previously opened.

% Revised September 2017 for new features in Vantage 3.3.0 SW release

%% initialize, clear variables, load .mat file

if ~simulateOnly
    % if "simulateOnly" is true, the system SW will not attempt to interact
    % with the HW at all, and operation will be limited to simulation mode
    % only.
    result = hardwareClose(); % close the hardware, in case it was left open from previous run.
    if ~strcmpi(result, 'Success')
        error(['VSX: unexpected result from hardwareClose: ', result]);
    end
    clear result
end

clear runAcq

% temporary cleanup of any spectral Doppler processing.
try
    spectralDoppler('cleanup');
catch
end

% Clear all variables with a few exceptions.
%
% NOTE:  "Mcr_" prefixed variables are designed for use by the MATLAB Compiler
% Runtime, though you don't have to use the MCR to use the "Mcr_" prefix.
vars = whos;
for i = 1:size(vars,1)
    if ~(strcmp(vars(i).name,'filename') || ...
         ~isempty(strfind(vars(i).name,'Mcr_')) || ... % Name contains "Mcr_"
         strcmp(vars(i).name,'RcvData') || ...
         strcmp(vars(i).name,'vars') || ...
         strcmp(vars(i).name,'Vs_VdasDebugging'))
        clear(vars(i).name);
    end
end

clear i vars



% call showEULA to determine if the EULA acceptance requirement has been
% met; a return value of 1 means it has
if showEULA == 0
    % user did not accept EULA so notify them and then quit here without doing anything else;
    fprintf(2, 'You must accept the terms of the EULA to use the Verasonics software.\n');
    clear
    return
end


% Read in structures from .mat file.
%
% Use filename variable if it is given and not empty, but only once.
% Otherwise ask the user for the filename.
if(~exist('filename', 'var') || (exist('filename', 'var') && isempty(filename)))
    % OK, the filename variable does NOT exist, so ask the user for the
    % filename.
    filename = input('Name of .mat file to process: ','s');
    if isempty(filename)
        disp('VSX exiting; no .mat file specified.')
        clear
        return
    end
end
% Use a try-catch so if file does not exist we can clear 'filename' for
% next run of VSX.
try
    load(filename);
    displayWindowTitle = filename;
    clear filename;
catch
    % File does not exist.  Clear 'filename' and report error.
    fprintf(2, ['The file "' filename '" was not found.\n']);
    clear filename;
    error('VSX: file not found');
end

%% set constants, defaults, & check required variables

% system constants
sysClk = 250; % system master oscillator clock rate in MHz

maxADRate = 62.5; % maximum A/D sample rate in MHz
% AFE5808 A/D maximum sample rate is 65 MHz.  So maximum for this system is
% 62.5 MHz (250 MHz sysClk divided by 4)
minADRate = 10; % Minimum A/D sample rate for AFE5808 in MHz

rchPerBd = 64; % number of receive channels per acquisition board
tchPerBd = 64; % number of transmit channels per acquisition board

%% Check for minimum required structures present, & set some defaults.
vars = whos;
RequiredStructs = {'Trans','Resource','TW','TX','Event'};
n = 0;
for i = 1:size(vars,1)
    if any(strcmp(vars(i).name, RequiredStructs)), n=n+1; end
end
if (n~=5), error('VSX: Trans, Resource, TW, TX and Event structures required as minimum. Exiting...\n'); end

% ***verbose level***:  set the variable to control command line display of
% informative or warning messages, if not set in user's script:
    % supported values:
    % verbose = 0: display error messages only
    % verbose = 1: display error and warning messages only
    % verbose = 2: display error, warning, and status messages
    % verbose = 3: display error, warning, status, and debug messages

if ~isfield(Resource, 'Parameters') || ~isfield(Resource.Parameters,'verbose') || isempty(Resource.Parameters.verbose)
    Resource.Parameters.verbose = 2; % default setting (see above definitions)
end


% numTransmit and numRcvChannels are required; report error if not found or
% out of range
% Simulate-only scripts can use up to 1024 channels but will crash matlab
% if more than 1024 are used, due to fixed memory allocation limits in
% runAcq.
if ~isfield(Resource.Parameters,'numTransmit') || ~isfield(Resource.Parameters,'numRcvChannels')...
        || isempty(Resource.Parameters.numTransmit) || isempty(Resource.Parameters.numRcvChannels)
    error('VSX: Resource.Parameters.numTransmit and/or Resource.Parameters.numRcvChannels not found.');
elseif Resource.Parameters.numTransmit < 1 || Resource.Parameters.numTransmit > 1024
    error('VSX: Resource.Parameters.numTransmit must be in the range 1 to 1024.');
elseif Resource.Parameters.numRcvChannels < 1 || Resource.Parameters.numRcvChannels > 1024
    error('VSX: Resource.Parameters.numRcvChannels must be in the range 1 to 1024.');
end


% - Set a Resource.Parameters.Connector default, if none provided, to
% identify which connector to use on HW configurations with more than one.
% For backward compatibility, check for 'connector' and if present translate it to
% 'Connector'.
if ~isfield(Resource.Parameters,'Connector') || isempty(Resource.Parameters.Connector)
    % check for old format before setting default
    if ~isfield(Resource.Parameters,'connector') || isempty(Resource.Parameters.connector)
        Resource.Parameters.Connector = 1; % default to connector #1
        % note this default value is implicit if system only has one connector,
        % but the variable still needs to be initialized since it is used in
        % system initialization logic for setting other variables.
        if Resource.Parameters.verbose > 2
            disp(' ');
            disp('VSX status: Resource.Parameters.Connector undefined; setting default value of 1.');
            disp(' ');
        end
    else
        % 'connector' exists so use it to define new format 'Connector',
        % but check for length = 1 to make sure it was not a 'Connector'
        % format definition with the name mis-spelled.
        if length(Resource.Parameters.connector) == 1 && Resource.Parameters.connector > 0
            % only one connector being used, so copy the value
            Resource.Parameters.Connector = Resource.Parameters.connector;
        else
            % connector = 0 so select both connectors (the old 'connector'
            % format was only used with UTA modules having one or two
            % connectors)
            Resource.Parameters.Connector = [1 2];
        end
        if Resource.Parameters.verbose > 2
            disp(' ');
            disp('VSX status: Resource.Parameters.Connector undefined; for backward compatibility');
            disp('    it is being created from older format Resource.Parameters.connector.');
            disp(' ');
        end
    end
end

% put entries in ascending order, to make decoding simpler
Resource.Parameters.Connector = sort(Resource.Parameters.Connector);

% for backward compatibility, don't allow deselecting all connectors
if isequal(Resource.Parameters.Connector, 0)
    Resource.Parameters.Connector = 1; % use default behavior of selecting connector 1
end

if ~isfield(Resource.Parameters, 'speedOfSound') || ...
        isempty(Resource.Parameters.speedOfSound)
    Resource.Parameters.speedOfSound = 1540;
end

if ~isfield(Resource.Parameters, 'speedCorrectionFactor') || ...
        isempty(Resource.Parameters.speedCorrectionFactor)
    Resource.Parameters.speedCorrectionFactor = 1.0;
end

if ~isfield(Resource.Parameters,'startEvent') || ...
        isempty(Resource.Parameters.startEvent)
    Resource.Parameters.startEvent = 1;
end

if ~isfield(Resource.Parameters, 'initializeOnly') || ...
        isempty(Resource.Parameters.initializeOnly)
    Resource.Parameters.initializeOnly = 0;
end

%% ***** Trans structure *****


% Check Trans structure for Trans.connType, or set default if not provided
if ~isfield(Trans, 'connType') || isempty(Trans.connType)
    % attempt to set a default of HDI-format connector(s) based on
    % Resource.Parameters.numTransmit (defaulting to HDI connectors
    % preserves backwards compatibility for pre-UTA scripts, since prior to
    % UTA the HDI connector was the only one available).
    if Resource.Parameters.numTransmit == 128 || (Resource.Parameters.numTransmit == 256 && isequal(Resource.Parameters.Connector, [1 2]))
        % HDI-format connector; note we must set this default value to
        % preserve backward compatibility with pre-UTA scripts, where the
        % HDI connector was the only one supported
        Trans.connType = 1;
        if Resource.Parameters.verbose>1
            disp(' ');
            disp('VSX status: Trans.connType not specified; setting a default value of 1 for use with HDI-format connectors.');
            disp(' ');
        end
    else
        % default to connType value of zero indicating simulation-only
        % operation with no UTA constraints applied.  If the user wants to
        % run with a HW system and a connector type other than 1, they must
        % specify the connType field in their setup script.  Note that for
        % test scripts intended to run on any HW configuration, connType
        % can be set to -1 meaning the script will automatically adapt to
        % whatever HW UTA module is present.
        Trans.connType = 0; % simulation-only operation with no UTA constraints.
        if Resource.Parameters.verbose>1
            disp(' ');
            disp('VSX status: Trans.connType not specified; setting a default value of 0');
            disp('for simulation-only operation with no UTA or system HW constraints applied.');
            disp(' ');
        end
    end
end

% Verify name is specified and set default if not
if ~isfield(Trans, 'name') || isempty(Trans.name)
    Trans.name = 'undefined'; % if user did not provide a name set to 'undefined' by default
    if Resource.Parameters.verbose > 1
        disp('VSX: Trans.name not specified.  Setting default name of ''undefined''.');
    end
elseif ~ischar(Trans.name)
    error('VSX: Trans.name must be a string variable.')
end

% Verify id is specified and set default if not
if ~isfield(Trans, 'id') || isempty(Trans.id)
    Trans.id = 0; % if user did not provide an id set to zero by default
    if Resource.Parameters.verbose > 1
        disp('VSX: Trans.id not specified.  Setting default value of zero.');
    end
elseif ischar(Trans.id)
    error('VSX: Trans.id must be a numeric id value, not a string variable.')
end

% Check for Trans.units and supply default if needed.  The default value of
% 'mm' is a compromise to V1 backward compatibility, but the warning
% message will alert the user to how to correct it.
if ~isfield(Trans, 'units')
    if Resource.Parameters.verbose
        fprintf(2,'VSX Warning:  Trans.units not specified.  Setting default units of mm.\n');
        fprintf(2,'If you are actually using wavelength units from the Trans structure \n');
        fprintf(2,'you must specify Trans.units = ''wavelengths'' in your script.\n');
    end
    Trans.units = 'mm';
end

% Verify units has a recognized value
if ~strcmp(Trans.units, 'wavelengths') && ~strcmp(Trans.units, 'mm')
    error('VSX: Unrecognized value for Trans.units. Must be ''mm'' or ''wavelengths''.')
end

% Verify frequency is specified and quit with an error if not
if ~isfield(Trans, 'frequency') || isempty(Trans.frequency)
    error('VSX: Trans.frequency must be specified in user script.')
end

% provide default bandwidth if not specified
if ~isfield(Trans, 'Bandwidth') || isempty(Trans.Bandwidth)
    if isfield(Trans, 'bandwidth') && ~isempty(Trans.bandwidth)
        % for backward compatibility, use the scalar bandwidth if it was
        % set by user in setup script
        Trans.Bandwidth = [-0.5, +0.5]*Trans.bandwidth + Trans.frequency;
        if Resource.Parameters.verbose > 1
            disp('VSX: Setting Trans.Bandwidth based on legacy Trans.bandwidth value.');
        end
    else
        Trans.Bandwidth = [0.7, 1.3] * Trans.frequency;  % 60% bandwidth default value
        if Resource.Parameters.verbose > 1
            disp('VSX: Trans.Bandwidth not specified. Setting default 60 percent value.');
        end
    end
end

% Verify numelements is specified and quit with an error if not
if ~isfield(Trans, 'numelements') || isempty(Trans.numelements)
    error('VSX: Trans.numelements must be specified in user script.')
end

% Verify ElementPos is specified and quit with an error if not
if ~isfield(Trans, 'ElementPos') || isempty(Trans.ElementPos)
    error('VSX: Trans.ElementPos must be specified in user script.')
end

% check for Trans.impedance and set default with warning if not provided
if ~isfield(Trans, 'impedance') || isempty(Trans.impedance)
    Trans.impedance = 20; % set default if not provided
    % The default value of 5 Ohms is absurdly low, and is intended to allow
    % the system to run but with an extremely low HV limit, thereby making
    % it easy for the user to discover the source of the restriction.
    if Resource.Parameters.verbose
        fprintf(2, 'VSX WARNING: A value of Trans.impedance has not been specified.  Default value of 20 Ohms has been set.\n');
        fprintf('This will allow script to run, but with extremely restricted transmit voltage limit.\n');
        fprintf('Actual Trans.impedance value should be added in user''s script.\n\n');
    end
end

% define default and max HV limits based on transducer
hvMax = computeTrans(Trans.name, 'maxHighVoltage');
hvDefault = 10*round(hvMax/20);

% - Check for maximum high voltage limit provided in Trans structure.  If not given, set a default.
if ~isfield(Trans,'maxHighVoltage') || isempty(Trans.maxHighVoltage)
    Trans.maxHighVoltage = hvDefault; % Transducer-specific default limit is defined above
    if Resource.Parameters.verbose>1
        fprintf(2,['VSX:  Trans.maxHighVoltage not specified.  Setting default limit of ', num2str(hvDefault), ' Volts.\n']);
    end
elseif Trans.maxHighVoltage < 1.6 || Trans.maxHighVoltage > hvMax
    % check if user-specified voltage is outside system HW or transducer maximum limits
    error(['VSX: Trans.maxHighVoltage must be within the range 1.6 to ', num2str(hvMax), ' Volts for ', Trans.name, ' transducer.'])
end

% Verify ElementSens is specified and quit with an error if not
if ~isfield(Trans, 'ElementSens') || isempty(Trans.ElementSens)
    error('VSX: Trans.ElementSens must be specified in user script.')
end

% Verify type is specified and quit with an error if not
if ~isfield(Trans, 'type') || isempty(Trans.type)
    error('VSX: Trans.type must be specified in user script.')
end

% Check for elBias field and set default if needed
if ~isfield(Trans, 'elBias') || isempty(Trans.elBias)
    % if not specified create the field and set it to zero (element bias
    % disabled)
    Trans.elBias = 0;
end

% Verify radius is specified for curved arrays and quit with an error if not
if Trans.type == 1 && (~isfield(Trans, 'radius') || isempty(Trans.radius))
    error('VSX: For curved arrays Trans.radius must be specified in user script.')
end

%% Check for Resource.Parameters.simulateMode and set default if needed
if ~isfield(Resource.Parameters,'simulateMode') || isempty(Resource.Parameters.simulateMode)
    if Trans.connType
        Resource.Parameters.simulateMode = 0; % default to running with HW if Trans.connType is not zero
    else
        Resource.Parameters.simulateMode = 1; % default to simulation-only if Trans.connType set to zero.
    end
    if Resource.Parameters.verbose > 1
        disp('VSX status: Resource.Parameters.simulateMode was not defined.');
        disp(['Setting a default value of ', num2str(Resource.Parameters.simulateMode), '.']);
    end
elseif Trans.connType == 0 && Resource.Parameters.simulateMode == 0
    % Trans.connType set to zero forces simulate-only operation (but note
    % simulateMode 2 is valid since it can be used with simulate-only
    % scripts)
    Resource.Parameters.simulateMode = 1;
    if Resource.Parameters.verbose > 2
        disp(' ');
        disp('VSX status: Resource.Parameters.simulateMode is being set to 1 for');
        disp('  simulation-only operation, since Trans.connType is set to zero.');
    end
end


% If .mat file contains an External Function definition, create the function in the temp directory
% and add the directory to the path.  Decoding the function here means that it can be used for a
% VsUpdate or GUI function.
if exist('EF','var')
    for i = 1:size(EF,2)
        if ~isfield(EF(i),'Function') || isempty(EF(i).Function), continue, end
        n = strfind(EF(i).Function{1},'=');
        if ~isempty(n)
            fname = textscan(EF(i).Function{1},'%*s %s %*[^\n]','Delimiter',{'=','('});
        else
            fname = textscan(EF(i).Function{1},'%s %*[^\n]', 'Delimiter','(');
        end
        fid = fopen([tempdir,fname{1}{1},'.m'], 'w');
        fprintf(fid,'function %s\n', EF(i).Function{1});
        for j = 2:size(EF(i).Function,2)
            fprintf(fid,'%s\n', EF(i).Function{j});
        end
    end
    status = fclose(fid);
    addpath(tempdir);
    rehash path
end


%% Open HW, Find HW configuration, check license status
% call the hwConfigCheck function to determine if HW system is present, and
% if so to determine the actual system configuration and check it for
% validity:
if Resource.Parameters.verbose > 1
    hccMode = 1; % enable verbose messages from hwConfigCheck
else
    hccMode = 0; % disable verbose messages from hwConfigCheck
end
Resource.SysConfig = hwConfigCheck(hccMode); % verbose as defined above
clear hccMode

% Get presence of HW system from hwConfigCheck results.
VDAS = Resource.SysConfig.VDAS;
numBoards = nnz(Resource.SysConfig.AcqSlots); % number of boards actually present, regardless of where they are

% If any configuration faults were found, display a warning message but
% ignore all configuration faults if HW test flag is set in the base
% workspace
if evalin('base', '~exist(''hwTestFlag'', ''var'')')
    if Resource.SysConfig.SWconfigFault
        % SW configuration fault means we cannot trust any HW or FPGA faults that
        % were detected, so ignore them and don't allow HW operation
        Resource.SysConfig.HWconfigFault = [];
        Resource.SysConfig.FPGAconfigFault = [];
        if Resource.Parameters.verbose
            fprintf(2, 'VSX WARNING: System is not in a valid released SW configuration.\n');
            fprintf(2, 'This probably means SW installation was not completed successfully.\n');
            if VDAS
                % don't display this text if HW isn't present, making it
                % meaningless
                fprintf(2, 'In this condition, HW and FPGA configuration checks may yield incorrect\n');
                fprintf(2, 'results, so VSX will be restricted to SW simulation use only.\n');
                fprintf(2, 'System can still be used with HW-SW diagnostic test and debug utilities.\n');
            end
            fprintf(2, 'Contact Verasonics Customer Support for additional information.\n');
            disp(' ');
        end
        VDAS = 0; % force to zero if HW was present, to prevent use of the HW
    elseif Resource.SysConfig.HWconfigFault > 1
        % HW is present but we can't communicate with it so just run in
        % simulation only mode
        Resource.SysConfig.FPGAconfigFault = []; % ignore any FPGA status information
        VDAS = 0; % force simulation only
        if Resource.Parameters.verbose
            fprintf(2, 'VSX WARNING: A HW system has been detected, but system SW is unable\n');
            fprintf(2, 'to communicate with it.  VSX will be restricted to SW simulation use only.\n');
            if Resource.SysConfig.HWconfigFault == 2
                fprintf(2, 'HW system is in Recovery Mode.\n');
            elseif Resource.SysConfig.HWconfigFault == 3
                fprintf(2, 'HW system is in an unidentified fault condition.\n');
            elseif Resource.SysConfig.HWconfigFault == 4
                fprintf(2, 'UTA adapter module is not installed or is in a fault condition.\n');
            elseif Resource.SysConfig.HWconfigFault == 5
                fprintf(2, 'HW system overheated, or fans not functioning.\n');
            end
            fprintf(2, 'System can still be used with HW-SW diagnostic test and debug utilities.\n');
            fprintf(2, 'Contact Verasonics Customer Support for additional information.\n');
            disp(' ');
        end
    elseif Resource.SysConfig.HWconfigFault == 1
        % HW is present but in unrecognized configuration
        if Resource.Parameters.verbose
            fprintf(2, 'VSX WARNING: A HW system has been detected, but it is an unrecognized\n');
            fprintf(2, 'system HW configuration.  VSX will be allowed to continue, but\n');
            fprintf(2, 'user scripts may not function properly.\n');
            fprintf(2, 'System can still be used with HW-SW diagnostic test and debug utilities.\n');
            fprintf(2, 'Contact Verasonics Customer Support for additional information.\n');
            disp(' ');
        end
    end

    if Resource.SysConfig.FPGAconfigFault
        fprintf(2, 'VSX WARNING: An out-of-date or unrecognized version of FPGA code\n');
        fprintf(2, 'has been found in the HW system.\n');
        disp('Enter "F" or just a carriage return to check and reprogram');
        disp(' all FPGA flash memories with the released code, (which may then');
        disp(' require a full power shutdown cycle before the system can be used);');
        reply = input('or enter a "Q" to exit VSX without reprogramming: ', 's');
        if strcmpi(reply, 'Q')
            % Quit VSX without doing anything
            return
        elseif strcmp(reply, 'proceed')
            % allow VSX to continue with existing FPGA code
        else
            % any other response means start FPGA code reprogramming using
            % hwConfigCheck mode 2
            [~] = hwConfigCheck(2);
            % Note that if hwConfigCheck actually reprograms the CGD or
            % ASC, it will force an exit with an error message informing
            % the user they must shutdown with a power cycle to
            % reinitialize PCIE devices.  Otherwise, it will return
            % normally and we can proceed, but to ensure that the
            % reprogrammed versions match the expected, quit VSX to force a
            % re-run.
            disp('FPGA programming completed without the need to reboot, please re-run VSX.');
            return
        end
    end
end

%% Check Connector type required by script versus UTA configuration

% First, replace 'wild card' Trans.connType with actual UTA value if HW is
% present
if Trans.connType == -1
    if VDAS
        % HW is present, so set connType to match the actual UTA
        Trans.connType = Resource.SysConfig.UTAtype(2);
    else
        % no HW, so set connType to 1 (default to HDI connector for HW
        % simulation)
        Trans.connType = 1;
    end
end


% now that HW-SW configuration checks and license checks are complete,
% force VDAS back to zero for simulation-only scripts
if VDAS && (Trans.connType == 0 || Resource.Parameters.simulateMode == 1)
    VDAS = 0;
end

% - Determine state of the variables used to control system operating limits.
if VDAS
    % running with HW system, so check state of HW configuration and options
    % Vantage 64 LE and TPC profile 5: use 'req' flags from hw as-is (note
    % that the p5req value has already been updated by hwConfigCheck to
    % reflect license file status).
    v64ena = Resource.SysConfig.v64req;
    p5ena = Resource.SysConfig.p5req;
else
    % HW is not present so set simulation defaults.  Note we can always use
    % the default SysConfig values returned by hwConfigCheck when no HW is
    % present, of all four boards being present.  The simulation SW driven
    % by the computeUTA function for the UTA required by the script will
    % only enable those CG's that are actually being used- any CG's that
    % are present but unused will be disabled, just as is the case on the
    % HW system.
    if Resource.Parameters.numRcvChannels == 64 && Resource.Parameters.numTransmit == 128
        % This is a script written explicitly for Vantage 64 LE so set the
        % flag identifying that.
        v64ena = 1;
    else
        v64ena = 0;
    end
    % if HW is present but VDAS has been forced to zero (such as in the
    % case of a system configuration problem), the Resurce.SysConfig fields
    % may have been set to states incompatible with simulation, so
    % overwrite those here with the same simulation defaults hwConfigCheck
    % would have returned if no HW was present:
    Resource.SysConfig.TXindex = 1; % standard frequency range
    Resource.SysConfig.RXindex = 5; % production AFE5812 (applies to all frequency ranges)
    % Check for optional field Resource.SimulationDefaults; if it exists,
    % use values from it to set the HW state to be used for simulation
    if isfield(Resource, 'SimulationDefaults')
        if isfield(Resource.SimulationDefaults, 'AcqSlots')
            Resource.SysConfig.AcqSlots = Resource.SimulationDefaults.AcqSlots;
            numBoards = nnz(Resource.SysConfig.AcqSlots);
        end
        % check first for user specification of the high-level
        % acqModuleConfig
        if isfield(Resource.SimulationDefaults, 'acqModuleConfig')
            switch Resource.SimulationDefaults.acqModuleConfig
                case {'LF', 'Low Frequency', 'Low'}
                    Resource.SysConfig.TXindex = 3;
                    Resource.SysConfig.RXindex = 5; % production AFE5812
                case {'SF', 'Standard Frequency', 'Standard'}
                    Resource.SysConfig.TXindex = 1;
                    Resource.SysConfig.RXindex = 5; % production AFE5812
                case {'HF', 'High Frequency', 'High'}
                    Resource.SysConfig.TXindex = 4;
                    Resource.SysConfig.RXindex = 5; % production AFE5812
                case 'HIFU'
                    Resource.SysConfig.TXindex = 1;
                    Resource.SysConfig.RXindex = 5; % production AFE5812
                otherwise
                    error('VSX: unrecognized string value for Resource.SimulationDefaults.acqModuleConfig.')
            end
        end
        % now let an explicit command for TXindex override what we just
        % did
        if isfield(Resource.SimulationDefaults, 'TXindex')
            Resource.SysConfig.TXindex = Resource.SimulationDefaults.TXindex;
        end
        if isfield(Resource.SimulationDefaults, 'RXindex')
            Resource.SysConfig.RXindex = Resource.SimulationDefaults.RXindex;
        end
    end
end


% if script intends to use HW, check for UTA
% compatibility and define UTA-specific mapping arrays.
if VDAS
    % At this point we know script intends to use HW and HW is present, so
    % need to check whether UTA required by script matches UTA that is
    % actually present.
    if Resource.SysConfig.UTAtype(2) ~= Trans.connType
        % UTA does not match script requirements; this is an error
        % condition
        if Resource.Parameters.verbose
            fprintf(2, 'VSX warning: System UTA module does not match script requirements.\n');
            fprintf(2, 'Enter s to continue in simulation-only mode, or just enter return to exit:  ');
            r = input(' ', 's');
            if strcmpi(r, 's')
                Resource.Parameters.simulateMode = 1;
                VDAS = 0;
                % For simulate mode independent of HW, we need to
                % explicitly set appropriate UTA type based on connector
                % type as given by Trans.connType.  So do nothing here and
                % then the following if statement will find VDAS = 0 and
                % set the UTA type.
            else
                % user wants to quit
                return
            end
        else
            % verbose is zero; don't query for user input- just exit
            % with an error message
            error('VSX: System UTA module does not match script requirements.\n')
        end
    end
end
if VDAS == 0 && Trans.connType > 0
    % no HW present, but we have a non-zero connType so set UTAtype to
    % match the script.
    switch Trans.connType
        case 1 % HDI format 260-pin connector
            if any(Resource.Parameters.Connector == 2)
                % using connector 2 or both connectors, so we need dual
                % connector UTA
                Resource.SysConfig.UTAtype = [1, 1, 2, 0];
            else
                % only using one connector; need to decide if this is
                % Vantage 64 script requiring mux UTA or not
                if Resource.Parameters.numTransmit == 64 && isfield(Trans, 'HVMux')
                    % This is a Vantage 64 mux UTA script since the above
                    % two conditions are met
                    Resource.SysConfig.UTAtype = [1, 1, 1, 2];
                    % UTA 260-Mux used only with Vantage 64 configuration
                else
                    Resource.SysConfig.UTAtype = [1, 1, 1, 0];
                    % UTA 260-S for all other configurations
                end
            end
        case 2 % breakout board or custom adapter; only one UTA type exists
            Resource.SysConfig.UTAtype = [1, 2, 1, 0];
        case 3 % Cannon 360 pin ZIF Connector; only one UTA type exists
            Resource.SysConfig.UTAtype = [1, 3, 1, 0];
        case 4 % Verasonics 408 pin connector; only one UTA type exists
            % Note UTA type [1 4 4 0] is the adapter STE test fixture for
            % use with PVT, but PVT does not support simulation operation
            % so we don't allow for it here
            Resource.SysConfig.UTAtype = [1, 4, 1, 0];
        case 6 % Hypertac Connector; only one UTA type exists, the
            % UTA-160DH/32 Lemo supporting three connectors (Two
            % Hypertac plus array of 32 Lemo)
            Resource.SysConfig.UTAtype = [1, 6, 3, 1];
        case 7 % GE 408 pin connector; only one UTA type exists
            Resource.SysConfig.UTAtype = [1, 7, 1, 0];
        case 8 % 1024 Mux direct connect adapter; only one UTA type exists
            Resource.SysConfig.UTAtype = [1, 8, 1, 3];
        otherwise
            error('VSX: Unrecognized Trans.connType value of %d.', Trans.connType);
    end
    % Since we just overwrote the UTAtype, wipe out UTAname to avoid
    % confusion
    Resource.SysConfig.UTAname = 'Modified for simulation operation';
end

% Add TPC structures not specified by the user.
if ~exist('TPC','var'), TPC = []; end
for i = length(TPC)+1:5
    TPC(i).inUse = [];
end

% Find which profiles are in use.
if exist('SeqControl','var')
    for j=1:size(SeqControl,2)
        if strcmp(SeqControl(j).command,'setTPCProfile')
            TPC(SeqControl(j).argument).inUse = 1;
        end
    end
else
    SeqControl = [];
end
% Profile 1 is active by default, if none other active.
if ~any([TPC.inUse])
    TPC(1).inUse = 1;
end

% Set TPC structure default values.
minTpcVoltage = 1.6;
for i = 1:size(TPC,2)
    if TPC(i).inUse
        if ~isfield(TPC(i), 'hv') || isempty(TPC(i).hv)
            TPC(i).hv = minTpcVoltage; % Default system setting for hv at startup.
        end
        if ~isfield(TPC(i), 'maxHighVoltage') || isempty(TPC(i).maxHighVoltage)
            TPC(i).maxHighVoltage = Trans.maxHighVoltage;
        else
            if TPC(i).maxHighVoltage > Trans.maxHighVoltage
                if Resource.Parameters.verbose
                    fprintf(2,'VSX WARNING: TPC(%d).maxHighVoltage reduced to Trans.maxHighVoltage limit of %.1f volts.\n',i,Trans.maxHighVoltage);
                end
                TPC(i).maxHighVoltage = Trans.maxHighVoltage;
            end
        end
        % Always initialize highVoltageLimit to match maxHighVoltage.
        TPC(i).highVoltageLimit = TPC(i).maxHighVoltage;
    else
        TPC(i).inUse = 0;
        TPC(i).hv = minTpcVoltage;
        TPC(i).maxHighVoltage = minTpcVoltage;
        TPC(i).highVoltageLimit = minTpcVoltage;
    end
end
clear minTpcVoltage


%% Determine mode of operation.
% - The variable VDAS is used to indicate presence of hardware.
% - the variable VDASupdates is used to control matlab generation of VDAS variables.
% - When not in simulateMode = 1 (simulate acquisition data), the system
% will automatically switch to simulate if HW is not available:

% check simulateMode value
switch Resource.Parameters.simulateMode
    case 0 % user wants to use hardware.
        rloopButton = 0;  % used to set state of rcv data loop button on GUI.
        VDASupdates = 1; % Tell matlab functions to update VDAS variables, whether we have HW or not:
        if VDAS
            % HW system is actually present, so get ready to use it
            simButton = 0;    % used to set state of simulate button on GUI.
        else
            % HW is not available so revert to simulate mode
            if Resource.Parameters.verbose>1
                if simulateOnly
                    fprintf('SW configured for simulation use only; switching to simulate mode.\n')
                else
                    fprintf('HW system is not available; switching to simulate mode.\n')
                end
            end
            simButton = 1;    % used to set state of simulate button on GUI.
            Resource.Parameters.simulateMode = 1;
        end
    case 1 % user wants to simulate acquisition
        VDAS = 0;  % don't talk to HW even if it is present
        VDASupdates = 0; % don't generate the VDAS variables either
        % but DO generate VDAS variables if profile 5 is to be used
        if TPC(5).inUse
            VDASupdates = 1;
        end
        rloopButton = 0;
        simButton = 1;
    case 2 % user wants to run script with existing RcvData array.
        rloopButton = 1;
        if VDAS && Trans.connType
            % HW system is present and script wants to use it
            VDASupdates = 1;  % Generate VDAS variables for running with HW
            simButton = 0;
        else
            % no HW or simulation-only so run in simulation mode
            VDASupdates = 0;
            % but DO generate VDAS variables if profile 5 is to be used
            if TPC(5).inUse
                VDASupdates = 1;
            end
            VDAS = 0;
            simButton = 0;
        end
    otherwise
        error('VSX: unrecognized value for Resource.Parameters.simulateMode.');
end

%% Resource.VDAS initialize
% If VDASupdates=1, check for presence of Resource.VDAS attributes. If not found set defaults.
if VDASupdates
    if ~isfield(Resource,'VDAS')
        Resource.VDAS = struct();
    end
    if ~isfield(Resource.VDAS,'numTransmit')
        Resource.VDAS.numTransmit = tchPerBd * numBoards;
    end
    if ~isfield(Resource.VDAS,'numRcvChannels')
        Resource.VDAS.numRcvChannels = rchPerBd * numBoards;
    end
    if ~isfield(Resource.VDAS,'exportDelta')
        Resource.VDAS.exportDelta = 2048;
    end
    if ~isfield(Resource.VDAS,'testPattern')
        Resource.VDAS.testPattern = 0;
    end
    if ~isfield(Resource.VDAS,'testPatternDma')
        Resource.VDAS.testPatternDma = 0;
    end
    if ~isfield(Resource.VDAS,'halDebugLevel')
        Resource.VDAS.halDebugLevel = 0;
    end
    if ~isfield(Resource.VDAS,'dmaPrecompute')
        Resource.VDAS.dmaPrecompute = 1;
    end
    if ~isfield(Resource.VDAS,'dmaComputeFree')
        Resource.VDAS.dmaComputeFree = 0;
    end
    if ~isfield(Resource.VDAS,'dmaTimeout')
        Resource.VDAS.dmaTimeout = 1000;
    end
    if ~isfield(Resource.VDAS,'watchdogTimeout') || isempty(Resource.VDAS.watchdogTimeout)
        % watchdog timeout interval in milliseconds; allowed range 10:10000
        % (10 msec to 10 sec)
        Resource.VDAS.watchdogTimeout = 10000; % 10 second default value
    end
    if ~isfield(Resource.VDAS,'el2ChMapDisable')
        Resource.VDAS.el2ChMapDisable = 0;
    end
    if ~isfield(Resource.VDAS,'sysClk')
        Resource.VDAS.sysClk = sysClk;
    end
    if ~isfield(Resource.VDAS,'elBiasSel')
        Resource.VDAS.elBiasSel = 0; % default state is element bias disabled
    end
end

%% Create channel mapping arrays based on UTA type
% At this point we have already confirmed SysConfig.UTAtype is compatible
% with the script being run; now we check if it is compatible with the
% system HW configuration, and if so create the mapping arrays to memory
% columns and VDAS channels.  For each recognized UTA type, determine if
% HW configuration is compatible and exit with error if not.
% Set VDAS to indicate presence/ absence of HW (for simulation, V128 & V256
% will be changed if needed to match UTA configuration)

if Trans.connType
    % don't call compute UTA if Trans.connType is zero
    UTA = computeUTA(Resource.SysConfig.UTAtype, Resource.Parameters.Connector);
end
if ~VDASupdates
    % simulation-only script or an unrecognized UTA, so create required
    % size of receive buffer etc. based on Resource.Parameters.numTransmit
    % with fixed 1:1 mapping
    numRbufCols = Resource.Parameters.numTransmit;
else
    if isequal(Resource.SysConfig.UTAtype, [1 2 1 0]) || Resource.VDAS.el2ChMapDisable
        % special case of breakout board with no connector mapping; number
        % of active channels will be defined here to match physical number
        % of acquisition modules in the system.  This also applies when
        % Resource.VDAS.el2ChMapDisable is true; in this case the UTA numCh
        % and activeCG fields will be overwritten to enable all channels on
        % all boards
        UTA.numCh = tchPerBd * numBoards;
        UTA.activeCG = upsample(Resource.SysConfig.AcqSlots, 2) + ...
            upsample(Resource.SysConfig.AcqSlots, 2, 1);
    end
    % find the number of boards in use, and determine if compatible with HW
    % system configuration
    activeBoards = min((UTA.activeCG([1 3 5 7]) + UTA.activeCG([2 4 6 8])), 1);
    if min(Resource.SysConfig.AcqSlots - activeBoards) < 0
        error('VSX: UTA module and script being used are not compatible with this HW system.');
    end
    % build the cgEnaDma and Cch2Vch arrays
    Resource.VDAS.cgEnaDma = sum(UTA.activeCG .* [1 2 4 8 16 32 64 128]);
    Resource.VDAS.Cch2Vch = [];
    % build HWactiveCG based on boards actually in HW system:
    HWactiveCG = [];
    for i=1:4
        if Resource.SysConfig.AcqSlots(i)
            HWactiveCG = [HWactiveCG, UTA.activeCG((2*i-1):(2*i))];
        end
    end
    activeCGnums = find(HWactiveCG);
    for i=1:length(activeCGnums)
        Resource.VDAS.Cch2Vch = [Resource.VDAS.Cch2Vch, (32*(activeCGnums(i)-1) + (1:32))];
    end
    numRbufCols = UTA.numCh;
end


clear rchPerBd tchPerBd


%% Check Resource.Parameters for valid attributes and initialize

% Check for VsUpdate function handle definition, and create default if needed
if ~isfield(Resource.Parameters,'UpdateFunction')||isempty(Resource.Parameters.UpdateFunction)
    Resource.Parameters.UpdateFunction = 'VsUpdate';
end
updateh = str2func(Resource.Parameters.UpdateFunction); % create handle to update function.

if VDASupdates
    if ~isfield(Resource.Parameters,'ProbeConnectorLED') || ...
            isempty(Resource.Parameters.ProbeConnectorLED)
        defaultProbeConnectorLED = [1 1 1 1];
        if isfield(Resource.VDAS,'shiConnectorLights')
            if ~isempty(Resource.VDAS.shiConnectorLights)
                Resource.Parameters.ProbeConnectorLED = zeros(1, 4);
                for i=1:length(Resource.Parameters.ProbeConnectorLED)
                    Resource.Parameters.ProbeConnectorLED(i) = ...
                        bitget(Resource.VDAS.shiConnectorLights, i);
                end
            else
                Resource.Parameters.ProbeConnectorLED = defaultProbeConnectorLED;
            end
            % Issue warning for use of deprecated Resource.VDAS.shiConnectorLights
            warning(['Resource.VDAS.shiConnectorLights is deprecated.\n' ...
                     'Use Resource.Parameters.ProbeConnectorLED = ' ...
                     '[%d, %d, %d, %d]\n'], Resource.Parameters.ProbeConnectorLED(:))
        else
            Resource.Parameters.ProbeConnectorLED = defaultProbeConnectorLED;
        end
        clear defaultProbeConnectorLED
    end
    if ~isfield(Resource.Parameters,'ProbeThermistor') || ...
            isempty(Resource.Parameters.ProbeThermistor)
        Resource.Parameters.ProbeThermistor = ...
            repmat(struct('enable', 0, ...
                          'threshold', 0, ...
                          'reportOverThreshold', 0), 1, 2);
        if isfield(Resource.VDAS,'shiThermistors')
            if ~all(size(Resource.VDAS.shiThermistors) == [3, 2])
                error('Resource.VDAS.shiThermistors must be 3 x 2');
            end
            for i = 1:length(Resource.Parameters.ProbeThermistor)
                Resource.Parameters.ProbeThermistor(i).enable = ...
                    Resource.VDAS.shiThermistors(1, i);
                Resource.Parameters.ProbeThermistor(i).threshold = ...
                    Resource.VDAS.shiThermistors(2, i);
                Resource.Parameters.ProbeThermistor(i).reportOverThreshold = ...
                    Resource.VDAS.shiThermistors(3, i);
            end
            warning(['Resource.VDAS.shiThermistors is deprecated. ' ...
                     'Resource.Parameters.ProbeThermistor has been populated ' ...
                     'from Resource.VDAS.shiThermistors. Inspect the ' ...
                     'Resource.Parameters.ProbeThermistor structure and ' ...
                     'use those values in your script.'])
        end
    end
    if ~isfield(Resource.Parameters,'SystemLED') || ...
            isempty(Resource.Parameters.SystemLED)
        defaultSystemLED = {'running', 'paused', 'activeTxAndOrRx', 'transferToHostComplete'};
        if isfield(Resource.VDAS, 'shiLeds')
            if ~isempty(Resource.VDAS.shiLeds)
                Resource.Parameters.SystemLED = cell(1, 4);
                for i = 1:length(Resource.Parameters.SystemLED)
                    shiLed = bitand(bitshift(Resource.VDAS.shiLeds, -8*(i-1)), hex2dec('FF'));
                    switch shiLed
                        case hex2dec('00')
                            Resource.Parameters.SystemLED{i} = 'off';
                        case hex2dec('0D')
                            Resource.Parameters.SystemLED{i} = 'running';
                        case hex2dec('13')
                            Resource.Parameters.SystemLED{i} = 'starting';
                        case hex2dec('01')
                            Resource.Parameters.SystemLED{i} = 'paused';
                        case hex2dec('02')
                            Resource.Parameters.SystemLED{i} = 'activeTxAndOrRx';
                        case hex2dec('12')
                            Resource.Parameters.SystemLED{i} = 'activeTxAndOrRxProfile5';
                        case hex2dec('0B')
                            Resource.Parameters.SystemLED{i} = 'transferToHostComplete';
                        case hex2dec('03')
                            Resource.Parameters.SystemLED{i} = 'pausedOnTriggerIn';
                        case hex2dec('20')
                            Resource.Parameters.SystemLED{i} = 'pausedDmaWaitPrevious';
                        case hex2dec('21')
                            Resource.Parameters.SystemLED{i} = 'pausedDmaWaitForProcessing';
                        case hex2dec('11')
                            Resource.Parameters.SystemLED{i} = 'pausedOnMultiSysSync';
                        case hex2dec('22')
                            Resource.Parameters.SystemLED{i} = 'pausedSync';
                        case hex2dec('10')
                            Resource.Parameters.SystemLED{i} = 'externalBoxFault';
                        case hex2dec('04')
                            Resource.Parameters.SystemLED{i} = 'missedTimeToNextAcq';
                        case hex2dec('07')
                            Resource.Parameters.SystemLED{i} = 'missedTimeToNextEB';
                        otherwise
                            error('Unrecognized shiLed value %d.\n', shiLed);
                    end
                end
            else
                Resource.Parameters.SystemLED = defaultSystemLED;
            end
            % Issue warning for use of deprecated Resource.VDAS.shiLeds.
            warning(['Resource.VDAS.shiLeds is deprecated.\n', ...
                     'Use Resource.Parameters.SystemLED = ' ...
                     '[''%s'', ''%s'', ''%s'', ''%s'']\n'], ...
                    Resource.Parameters.SystemLED{:})
        else
            Resource.Parameters.SystemLED = defaultSystemLED;
        end
        clear defaultSystemLED shiLed
    end
end

% - If VDASupdates are requested, validate numTransmit and numRcvChannels.
%   When running with the hardware, Resource.Parameters.numTransmit and
%   Resource.Parameters.numRcvChannels must match the connector channel
%   count for the UTA module or SHI configuartion being used; for
%   simulation without VDASupdates, the match with hardware is not
%   enforced.
%
%   The number of columns in the Receive data buffer must also be
%   constrained to match the HW configuration, but only when running with
%   the HW.
if VDASupdates
    if (Resource.Parameters.numTransmit~=numRbufCols)
        error(['VSX: Resource.Parameters.numTransmit must equal ', num2str(numRbufCols), ...
            ', the number of connector channels for connector configuration being used.']);
    end
    if isfield(Resource, 'RcvBuffer')
        % verify column size in receive buffer matches numRbufCols, as
        % set above for the system configuration and UTA being used
        if Resource.RcvBuffer(1).colsPerFrame ~= numRbufCols
            % Need to modify RcvBuffer size & warn user
            if Resource.Parameters.verbose>1
                disp(' ');
                disp(['VSX Status: Resource.RcvBuffer.colsPerFrame is being changed to ', num2str(numRbufCols)]);
                disp(' to allow the script to run on the Vantage system HW configuration.');
                disp(' ');
            end
        end
        % There must be either 1 or an even number of frames in RcvBuffer
        for i=1:size(Resource.RcvBuffer,2)
            if Resource.RcvBuffer(i).numFrames ~= 1 && ...
                    mod(Resource.RcvBuffer(i).numFrames, 2)
                error(sprintf(['VSX: Resource.RcvBuffer(%d).numFrames ' ...
                               'must be 1 or an even integer.'], i));
            end
        end
    end
else
    % for simulation without VDASupdates, still need to enforce
    % Resource.RcvBuffer.colsPerFrame = numTransmit.
    if Resource.RcvBuffer(1).colsPerFrame ~= Resource.Parameters.numTransmit
        % Need to modify RcvBuffer size & warn user
        if Resource.Parameters.verbose>1
            disp(' ');
            disp(['VSX Status: Resource.RcvBuffer.colsPerFrame is being changed to ', num2str(Resource.Parameters.numTransmit)]);
            disp(' to allow the script to run with Vantage system SW.');
            disp(' ');
        end
    end
end
if isfield(Resource, 'RcvBuffer')
    % Now make column size correction to all RcvBuffers, but don't need to
    % repeat warnings. Also determine total defined RcvBuffer size and
    % report to user at verbose level 3
    rBufSize = 0;
    for i=1:size(Resource.RcvBuffer,2);
        Resource.RcvBuffer(i).colsPerFrame = Resource.Parameters.numTransmit;
        rBufSize = rBufSize + Resource.RcvBuffer(i).colsPerFrame * ...
            Resource.RcvBuffer(i).rowsPerFrame * Resource.RcvBuffer(i).numFrames;
    end
    if Resource.Parameters.verbose > 2
        disp(['VSX Debug: Total allocated Receive Buffer memory space is ' num2str(rBufSize/2^19) ' Megabytes.']);
    end
end

% - Check for the number of LogData records specified.  These records are use for debugging
%   purposes and are optionally generated by the mex file runAcq.  Default number is 128 records.
%   Meaning of fields in a 4 uint LogData record:
%     int32    id        This number is an identifier for the record and is
%                         the same as 'id' in LOGTIME(id).
%     int32    datatype  Identifier for type of data (0 = time, 1 = int, 2 = double);
%     int32     data1     For a time record, this field is the hi value of
%                         an AbsoluteTime number; for data record, it is
%                         the an integer value (for datatype=1) or the integer
%                         representation of a double (hi 16 bits integer,
%                         low 16 bits fraction) (for datatype=2).
%     int32     data2     For a time record, this field is the lo value of
%                         an AbsoluteTime number.

global LogData;
shouldConvertLogData = true;
if ~isfield(Resource.Parameters,'numLogDataRecs') || ...
        isempty(Resource.Parameters.numLogDataRecs)
    Resource.Parameters.numLogDataRecs = 128;
    shouldConvertLogData = false;
end
LogData = zeros(4, Resource.Parameters.numLogDataRecs, 'int32');




%% **** Determine whether to reuse RcvData buffer ****.
% If RcvData already exists in workspace, check for same dimensions
% as specified in Resource structure. If dimensions not equal, clear the RcvData.
if exist('RcvData','var')
    clrRcvData = 0;  % set default to not clear RcvData
    if size(RcvData,1) ~= size(Resource.RcvBuffer,2)
        clrRcvData = 1;
    else
        for i = 1:size(RcvData,1)
            if size(RcvData{i},1) ~= Resource.RcvBuffer(i).rowsPerFrame
                clrRcvData = 1;
            end
            if size(RcvData{i},2) ~= Resource.RcvBuffer(i).colsPerFrame
                clrRcvData = 1;
            end
            if size(RcvData{i},3) ~= Resource.RcvBuffer(i).numFrames
                clrRcvData = 1;
            end
            if isfield(Resource.RcvBuffer,'datatype')
                if ~strcmp(class(RcvData{i}),Resource.RcvBuffer(i).datatype)
                    fprintf('RcvData found in workplace, but datatype is not ''int16'' - clearing.\n');
                    clrRcvData = 1;
                end
            elseif ~strcmp(class(RcvData{i}),'int16')
                fprintf('RcvData found in workplace, but datatype is not ''int16'' - clearing.\n');
                clrRcvData = 1;
            end
        end
    end
    if (clrRcvData == 0)
        fprintf('RcvData in workplace matches Resource.RcvBuffer specification - reusing without clearing.\n');
    else
        if Resource.Parameters.simulateMode==2
            fprintf(2,'SimulateMode 2 specified, but RcvData in workspace doesn''t match definition in Resource.RcvBuffer.\n');
        end
        clear('RcvData');
    end
end


%% ***** Trans.HVMux and .Connector mapping *****

% Check for the special case of the Vantage 64 HV Mux adapter
if isequal(Resource.SysConfig.UTAtype, [1, 1, 1, 2])
    % The mux adapter only works with non-mux scripts
    if computeTrans(Trans.name,'HVMux')
        error('VSX: Probe with HVMux cannot be used with UTA 260-MUX adapter module.');
    end
    if ~isfield(Trans,'HVMux') || isempty(Trans.HVMux)
        error('VSX: Trans.HVMux structure required for UTA 260-MUX.  Use computeUTAMux64 in your script.');
    end
end


if isfield(Trans,'HVMux')
    % This is an HVMux transducer
    % check for settling time field and create with default value if not
    % present
    if ~isfield(Trans.HVMux,'settlingTime') || isempty(Trans.HVMux.settlingTime)
        Trans.HVMux.settlingTime = 4; % default to 4 usec
    end
    % verify columns are correct length in Trans.HVMux.Aperture
    if size(Trans.HVMux.Aperture, 1) ~= Trans.numelements
        error('VSX: Number of rows in Trans.HVMux.Aperture must equal Trans.numelements.')
    end
    sizeTransActive = nnz(Trans.HVMux.Aperture(:,1)); % use first column to set sizeTransActive
    for i=1:size(Trans.HVMux.Aperture, 2) % now check all other columns for equivalent value
        if sizeTransActive ~= nnz(Trans.HVMux.Aperture(:,i))
            error('VSX: All columns in Trans.HVMux.Aperture must have the same size active aperture.')
        end
    end
    % now check for illegal channel index in the Aperture array
    if max(Trans.HVMux.Aperture(:)) > Resource.Parameters.numTransmit
        error('VSX: Trans.HVMux.Aperture is indexing a non-existent system channel number.\n')
    end
elseif isfield(Trans,'Connector')
    % No HVMux so we use Trans.Connector instead
    if size(Trans.Connector, 1) ~= Trans.numelements
        % verify column of correct length
        error('VSX: Length of Trans.Connector must equal Trans.numelements.')
    end
    % now check for illegal channel index in the Aperture array
    if max(Trans.Connector(:)) > Resource.Parameters.numTransmit
        error('VSX: Trans.Connector is indexing a non-existent system channel number.\n')
    end
    sizeTransActive = nnz(Trans.Connector); % use first column to set sizeTransActive
else
    % No Trans.Connector either so create a default version with 1:1 mapping if possible
    if Trans.numelements ~= Resource.Parameters.numTransmit
        error('VSX: Cannot create default Trans.Connector array since Trans.numelements ~= Resource.Parameters.numTransmit.')
    end
    Trans.Connector = (1:Trans.numelements)';
    sizeTransActive = Trans.numelements;
end

if VDASupdates && Resource.VDAS.el2ChMapDisable
    % when this flag is set while running with HW, overwrite
    % Trans.numelements and Trans.Connector with the forced values of all
    % HW channels and a one-to-one mapping
    if ~isfield(Trans,'HVMux')
        Trans.numelements = UTA.numCh;
        sizeTransActive = Trans.numelements;
        if isfield(Trans,'Connector')
            Trans.Connector = (1:Trans.numelements)';
        end
    end
elseif isequal(Resource.SysConfig.UTAtype, [1, 6, 3, 1])
    % for the UTA 160-DH/32 Lemo we must now "re-map" the Trans.Connector
    % array with the UTA mapping as defined by the computeUTA function.
    Trans.Connector = UTA.TransConnector(Trans.Connector);
end

%% ***** Media structure *****
% - Check for existance of 'MP' structure array. If not found, set defaults.
%       Media points are specified by number, position(x,y,z), and reflectivity.
%       MP(1,:) = [0,0,20,1.0];
if ~exist('Media','var')
    Media.model = 'PointTargets1';
    pt1;
    Media.numPoints = size(Media.MP,1);
elseif (isfield(Media, 'model'))
    if (strcmp(Media.model, 'PointTargets1'))
        pt1;
    elseif (strcmp(Media.model, 'PointTargets2'))
        pt2;
    elseif (strcmp(Media.model, 'PointTargets3'))
        mpt3;
    else
        error('Unknown media model. Could not initialize MP array.\n');
    end
elseif (isfield(Media, 'program'))
    eval(Media.program);
    Media.numPoints = size(Media.MP,1);
elseif (~isfield(Media, 'MP'))
    error('MP array not found.\n');
elseif (~isfield(Media, 'numPoints'))
    Media.numPoints = size(Media.MP,1);
end


%% ***** PData structure *****
if exist('PData','var')
    % - Check to see if the PData structure(s) specify all pdeltas and Regions.  If not, create them.
    for i = 1:size(PData,2)
        k = 0;  % k is flag for need to call computeRegions
        if (~isfield(PData(i),'Region'))||(isempty(PData(i).Region)) % if no Regions defined
            k = 1;  % computeRegions will create and compute a single Region the size of PData
        else  % some Regions are defined, but may not be computed.
            for j = 1:size(PData(i).Region,2)  % check for all Regions computed.
                if (~isfield(PData(i).Region(j),'numPixels'))||(~isfield(PData(i).Region(j),'PixelsLA'))
                    % no compute fields defined - need to specify Shape structure.
                    if (~isfield(PData(i).Region(j),'Shape'))||(isempty(PData(i).Region(j).Shape))
                        error('No Shape structure specified in PData(%d).Region(%d)\n',i,j);
                    end
                    k = 1;  % both computed fields not found but Shape structure provided.
                elseif (isempty(PData(i).Region(j).numPixels))||(isempty(PData(i).Region(j).PixelsLA))
                    % computed fields found but one or both empty; in this case also need Shape attribute
                    if (~isfield(PData(i).Region(j),'Shape'))||(isempty(PData(i).Region(j).Shape))
                        error('No Shape structure specified in PData(%d).Region(%d)\n',i,j);
                    end
                    k = 1;
                end
            end
        end
        if (k==1), [PData(i).Region] = computeRegions(PData(i)); end  % compute all Regions
    end
    % - Capture original PData(1).pdeltas for zoom function.
    if isfield(PData(1),'PDelta')
        if size(PData(1).PDelta,2) ~= 3
            error('VSX: PData(1).PDelta array must have 3 values - X,Y and Z.');
        end
        orgPdeltaX = PData(1).PDelta(1);
        orgPdeltaZ = PData(1).PDelta(3);
    elseif isfield(PData(1),'pdelta')
        orgPdeltaX = PData(1).pdelta;
        orgPdeltaZ = PData(1).pdelta;
    else
        orgPdeltaX = PData(1).pdeltaX;
        orgPdeltaZ = PData(1).pdeltaZ;
    end
end


%% ***** HIFU Option Check *****
% Check for use of HIFU transmit profile 5 in setup script and initialize the workspace variable,
% 'TPC(5).inUse', which controls all HIFU or profile 5-related features in the system:
%     inUse = 0 (default) if the setup script does not make any use of Profile 5;
%     inUse = 1 for Extended Burst Option using internal auxiliary power supply for Profile 5 transmit.
%     inUse = 2 for HIFU semi-custom option using ext. power supply with remote control for Profile 5 transmit.
% Check SeqControl commands for a setTPCProfile command with an argument of 5.

% check for a user-specified substitute for TXEventCheck function & create
% handle, even if profile 5 is not in use at the moment (it may be after
% some gui activity, etc.)
if ~isfield(Resource,'HIFU') || ~isfield(Resource.HIFU,'TXEventCheckFunction') || isempty(Resource.HIFU.TXEventCheckFunction)
    Resource.HIFU.TXEventCheckFunction = 'TXEventCheck';
end
% now confirm the specified function actually exists and is on the path
if isempty(which(Resource.HIFU.TXEventCheckFunction))
    error('VSX: The function specified by ''Resource.HIFU.TXEventCHeckFunction'' could not be found.');
end
TXEventCheckh = str2func(Resource.HIFU.TXEventCheckFunction); % create handle to TXEventCheck function.



% If profile 5 is going to be used, check that the hardware has the appropriate HIFU or Extended Transmit option installed.
% Check for the Resource.Parameter attribute, 'externalHifuPwr' from the
% setup script indicating that the external power supply is to be used.
% Also check and initialize other parameters related to use of profile 5.
if (TPC(5).inUse == 1)

    % Now check for the required TPC(5) max high voltage limit
    if ~exist('TPC','var')||size(TPC,2)<5||~(isfield(TPC(5),'maxHighVoltage'))||isempty(TPC(5).maxHighVoltage)
        error('VSX: A script using profile 5 transmit must specify TPC(5).maxHighVoltage.');
    end

    if isfield(Trans, 'HVMux')
        % if a probe with HV mux chips is being used, exit with an error
        % unless user has set special field to enable use
        if isfield(Trans.HVMux, 'P5allowed') && Trans.HVMux.P5allowed == 1
            if Resource.Parameters.verbose
                fprintf(2, 'VSX WARNING: A Sequence Control Command selecting TPC Profile 5 has been detected.\n');
                fprintf(2, '             Use of extended burst durations can be destructive to an HVMux probe.\n');
            end
            if Resource.Parameters.verbose > 1
                ans = input('Enter "y" to contimue, or return to exit: ', 's');
                if ~strcmpi('y', ans)
                    return
                end
                clear ans
            end
        else
            error('VSX: Profile 5 transmit is not allowed when using probes with HVMux element switching.');
        end
    end

    % Now check actual HW configuration
    if VDAS
        % HW is present, so query actual HW configuration
        % Check for the 'p5ena flag. The value might be:
            % 0 for no HIFU option enabled
            % 1 for the "Extended Transmit Option" using internal auxilliary supply on the TPC
            % 2 for the "HIFU Option" using external OEM 1200 Watt power supply.
        if p5ena == 0
            error('VSX: Profile 5 transmit features can not be used since neither Extended Transmit or HIFU option is installed in the system.');
        end

        if p5ena == 2
            if isfield(Resource.HIFU,'externalHifuPwr') && (Resource.HIFU.externalHifuPwr >= 1)
                TPC(5).inUse = 2;
            else
                error(['VSX: System is configured for use with external HIFU power supply, ',...
                'but Resource.HIFU.externalHifuPwr has not been set in SetUp script.']);
            end

            if Resource.SysConfig.HWconfigFault == 1
                % HIFU operation not allowed if configuration faults are
                % present
                error('VSX: HIFU feature cannot be used due to HW configuration fault or licensing error.');
            end
            % At this point, the system is set to use the external power
            % supply option.  So we have to initialize the external power
            % supply and make sure it is connected and up to the initial
            % TPC(5).hv Voltage setting before we try to open the VDAS HW
            % (TPC initialization is likely to fail if HIFU capacitor
            % voltage is sitting at zero). - Initialize external supply
            % communication and set it to 1.6 Volts (or user-specified
            % TPC(5).hv value), 2 Amps

            % As part of initializing the power supply, we also determine
            % whether it is configured for series or parallel connection of
            % the two outputs and configure the control function to match.
            % This also lets us determine whether the power supply is
            % actually connected to the system and working properly, since
            % we will read back the actual voltage from the TPC.

            % first determine if push capacitor is still charged up from a
            % previous run of the system
            [Result,pushCapVoltage] = getHardwareProperty('TpcExtCapVoltage');
            if ~strcmp(Result,'Success')
                error('VSX: Error from getHardwareProperty call to read push capacitor Voltage.');
            end
            if pushCapVoltage > 2
                % push cap not fully discharged, so wait longer after
                % enabling external supply to let it bleed down to
                % programmed value
                psPause = 2; % wait 2 seconds
            else
                psPause = 0.7; % shorter time is enough otherwise
            end

            Resource.HIFU.extPwrConnection = 'parallel'; % assume it is parallel to start
            % now try programming the supply
            [extPSstatus, ~] = extPwrCtrl('INIT', TPC(5).hv);
            if extPSstatus ~= 0
                error('VSX: Cannot initialize external power supply.  Make sure it is connected and turned on.');
            end
            pause(psPause); % wait for the power supply to come up
            % now read the actual voltage at push cap. through the TPC
            [Result,pushCapVoltage] = getHardwareProperty('TpcExtCapVoltage');
            if ~strcmp(Result,'Success')
                error('VSX: Error from getHardwareProperty call to read push capacitor Voltage.');
            end

            % now check the result:
            if pushCapVoltage < 0.7 * TPC(5).hv
                % less than this Voltage means we have a fault or power supply
                % not connected;
                error('VSX: No output from external power supply, or it may be disconnected.');
            elseif pushCapVoltage >= 0.7 * TPC(5).hv && pushCapVoltage < 1.4 * TPC(5).hv
                % between these levels means supply is using parallel
                % connection and is working properly, so no need to do
                % anything else
            elseif pushCapVoltage >= 1.4 * TPC(5).hv && pushCapVoltage < 2.8 * TPC(5).hv
                % between these levels means we have a
                % series-connected supply, so reconfigure for series
                % operation and test again

                % but first we must disable supply so it can be
                % reprogrammed
                [~, ~] = extPwrCtrl('CLOSE',1.6);
                pause(0.5) % give it time to close
                Resource.HIFU.extPwrConnection = 'series'; % reconfigure for series
                % now try reprogramming the supply
                [extPSstatus, ~] = extPwrCtrl('INIT',TPC(5).hv);
                if extPSstatus ~= 0
                    error('VSX: Cannot initialize external power supply.  Make sure it is connected and turned on.');
                end
                pause(1.5); % wait 1.5 seconds for the power supply to come up & bleed capacitor down
                % now read the actual voltage at push cap. through the TPC
                [Result,pushCapVoltage] = getHardwareProperty('TpcExtCapVoltage');
                if ~strcmp(Result,'Success')
                    error('VSX: Error from getHardwareProperty call to read push capacitor Voltage.');
                end
                if pushCapVoltage < 0.7 * TPC(5).hv || pushCapVoltage > 1.4 * TPC(5).hv
                    error('VSX: Error while initializing external supply for series connection.');
                end
            else
                % any other voltage means we have a fault.
                error('VSX: Error while initializing external supply for parallel connection.');
            end

        else % system is configured for internal aux. supply; confirm that's what the user intended
            if isfield(Resource.HIFU,'externalHifuPwr') && (Resource.HIFU.externalHifuPwr >= 1)
                error('VSX: System is not configured for use with external HIFU power supply.');
            end
        end
    else
        % No hardware is present, so set HIFU defaults for simulation.
        if isfield(Resource.HIFU,'externalHifuPwr') && (Resource.HIFU.externalHifuPwr >= 1)
            % set simulation defaults for external HIFU
            TPC(5).inUse = 2;
        end
    end
end


%% ***** TW Structure *****

% For the TW structure all checks for required veriables, generation of
% default values for optional variables, and checking for required range
% limits and format of all variables are executed in the computeTWWaveform
% function which is called automatically by VsUpdate(TW).  The VsUpdate
% function is called by runAcq during initialization for the TW, TX, and
% Receive structures; VsUpdate is also called during run time from within runAcq,
% whenever an 'update&Run' command for TW or TX is executed.

%% ***** TX Structure *****

% Call the VsUpdate(TX) function to check for missing variables, out of range
% values, etc.  It will also add the TX VDAS parameters for running
% with the hardware, if VDASupdates is true.

if ~any(isfield(TX,{'VDASApod'}))
    % VsUpdate(TX) call is made by runAcq during initialization, so at this
    % point we do nothing and just continue
else
    fprintf(2,'VSX: TX.VDASApod cannot be included in SetUp file. To set TX.VDASApod\n');
    fprintf(2,'manually, use Resource.Parameters.updateFunction to specify user provided function.\n');
    error('VSX exiting...');
end



%% ***** DMAControl Structure *****

if exist('DMAControl', 'var')
    error(['%s: DMAControl structures cannot be included in SetUp file.\n', ...
           'To specify DMAControl structures, define an update functiion ', ...
           'using Resource.Parameters.updateFunction.\n'], mfilename)
end


%% **** RcvProfile Structure ***

if VDASupdates % only add structure when preparing to run with HW
    % if RcvProfile structure doesn't exist, create it before calling compute
    % function
    if ~exist('RcvProfile','var')
        RcvProfile.AntiAliasCutoff = [];
    end

    RcvProfile = computeRcvProfile(RcvProfile);
    % The compute function will assign defaults for unspecified items, check
    % for valid values, etc.
end


%% **** TGC Structure *****
% If a TGC(1) waveform exists, initialize the tgc slide pot variables.
if exist('TGC','var')
    tgc1 = double(TGC(1).CntrlPts(1));
    tgc2 = double(TGC(1).CntrlPts(2));
    tgc3 = double(TGC(1).CntrlPts(3));
    tgc4 = double(TGC(1).CntrlPts(4));
    tgc5 = double(TGC(1).CntrlPts(5));
    tgc6 = double(TGC(1).CntrlPts(6));
    tgc7 = double(TGC(1).CntrlPts(7));
    tgc8 = double(TGC(1).CntrlPts(8));
else
    tgc1=511; tgc2=511; tgc3=511; tgc4=511; tgc5=511; tgc6=511; tgc7=511; tgc8=511;
end


%% **** ReconInfo structures ****

% The functionality to automatically generate ReconInfo.Aperture is
% included in VsUpdate(Receive), which has already been called so nothing
% needs to be done here.  By locating this code in VsUpdate(Receive) it will
% always be executed even when the Receive structure is modified during
% system operation.


%% **** Recon Structure ****
% - Check 'Recon' structure array for compatible PData and destination buffer sizes, and number
%   of ReconInfos.
%   If only one column and ReconInfo.regionnum = 0 or [], add ReconInfo structures for all regions.
if exist('Recon','var')
    nextRI = size(ReconInfo,2) + 1; % this will keep track of new ReconInfo structures created.
    for i = 1:size(Recon,2)
        j = Recon(i).pdatanum;
        if isfield(Recon(i),'IntBufDest')&&(~isempty(Recon(i).IntBufDest))&&(Recon(i).IntBufDest(1) > 0)
            k = Recon(i).IntBufDest(1); % get dest. buffer no.
            % Check for InterBuffer specified.
            if ~isfield(Resource,'InterBuffer')
                error('Resource.InterBuffer(%d) specified in Recon but not defined.\n',k);
            end
        else
            % If no IntBufDest, check to verify that no IQ reconstructions are specified.
            for k = 1:numel(Recon(i).RINums)
                n = Recon(i).RINums(k);
                if (isfield(ReconInfo(n),'mode'))&&(~isempty(ReconInfo(n).mode))
                    if isnumeric(ReconInfo(n).mode)
                        if ReconInfo(n).mode > 2
                            error('VSX: Recon(%d) specifies a ReconInfo.mode > 2, but no IntBufDest is given.\n',i);
                        end
                    elseif ischar(ReconInfo(n).mode)
                        switch ReconInfo(n).mode
                            case {'replaceIntensity','addIntensity','multiplyIntensity'}
                                continue
                            otherwise
                                error('VSX: Recon(%d) specifies a ReconInfo.mode > 2, but no IntBufDest is given.\n',i);
                        end
                    else
                        error('VSX: Recon(%d).mode must be numeric or string.\n',i);
                    end
                else
                    error('VSX: Recon(%d).mode missing or empty.\n',i);
                end
            end
        end
        if isfield(Recon(i),'ImgBufDest')&&(~isempty(Recon(i).ImgBufDest))
            k = Recon(i).ImgBufDest(1); % get dest. buffer no.
            if k ~= 0
                % Check for ImageBuffer specified.
                if ~isfield(Resource,'ImageBuffer')
                    error('Resource.ImageBuffer(%d) specified in Recon but not defined.(%d)\n',k,j);
                end
            end
        end
        % Set the required value for Recon(i).numchannels.  Note that even
        % if a value was set by the user's script it will be overwritten
        % here.  This is to preserve backward compatibility to V1 scripts.
        % With new HAL CG mapping, Recon numchannels should always match
        % the size of the receive buffer:
        Recon(i).numchannels = Resource.RcvBuffer.colsPerFrame;
        % check for only one column and all regionnums in specified ReconInfos set to zero.
        if size(Recon(i).RINums,2)==1
            for j = 1:size(Recon(i).RINums,1) % set any missing or empty fields to zero.
                if ~isfield(ReconInfo(Recon(i).RINums(j,1)),'regionnum')||...
                        isempty(ReconInfo(Recon(i).RINums(j,1)).regionnum)
                    ReconInfo(Recon(i).RINums(j,1)).regionnum = 0;
                end
            end
            % check for all ReconInfos in col. set to specify all regions (0)
            if ~any([ReconInfo(Recon(i).RINums(:,1)).regionnum])
                for j = 1:size(Recon(i).RINums,1) % set first col RIs to region 1
                    ReconInfo(Recon(i).RINums(j,1)).regionnum = 1;
                end
            else
                if ~all([ReconInfo(Recon(i).RINums(:,1)).regionnum])
                    error('VSX: single col. Recon(%d).RINums must specify ReconInfo.regionnums that are all set or all missing or 0.\n',i)
                    return
                end
            end
        end
    end
end


%% If not in simulate mode, try to open the hardware.
if VDAS
    % Try to open the Verasonics hardware.
    %
    % NOTE:  Textual output from compiled C code is not seen in MATLAB until
    % after the function returns.
    try
        Result = hardwareOpen();
    catch
        Result = 'FAIL';
    end
    if strcmpi(Result,'SUCCESS')
        % Hardware opened successfully.
        %   When using a debugger, inform the Verasonic HAL so that it can
        %   disable Watch Dog timout, and maximize timeouts, such as for DMA.
        if exist('Vs_VdasDebugging', 'var')
            % OK, the flag exists, so set state accordingly.
            if(1 == Vs_VdasDebugging)
              % Warn user that Watch Dog timer is disabled.
              if Resource.Parameters.verbose
                  display('WARNING:  Vs_VdasDebugging set to 1, so DISABLING the VDAS watch dog timer and informing the HAL to be in debug mode.');
              end
            end
            Result = vdasDebugging(Vs_VdasDebugging);
            if ~strcmpi(Result, 'SUCCESS')
              % ERROR!  Failed to set debugging state.
              error('ERROR!  Failed to set Verasonics VDAS debugging state.');
            end
        end
    else  % If hardware didn't open successfully, report an error and quit
        error('VSX:  HW has been detected, but the hardwareOpen function failed.')
    end
end


%% check probe connected status and ID value
if VDAS == 1 && Resource.SysConfig.UTA
    % skip the probe ID and connector selection if there is no HW system
    % present, or if there is HW but no SHI or UTA baseboard is present
    % (which will be indicated by SysConfig.UTA set to zero)
    % Query the system for scanhead connectivity and select scanhead.

    % as a workaround for the VTS-321 problem where the first read after
    % power-up of the L22-8v returns an incorrect ID of zero, we do the probe
    % ID query twice and ignore the results the first time.
    [p1, s1, Id1] = vdasScanheadQuery();
    clear p1 s1 ID1
    [Presence, Selected, ID] = vdasScanheadQuery();
%     if length(Resource.Parameters.Connector) > 1
%         error('VSX: code for selecting more than one connector has not been updated yet.');
% %         if isequal(Resource.SysConfig.UTAtype, [1 1 2 0])
% %             % UTA 260-D is present
% %             % for special case of connector set to zero we use both connectors
% %             % to support 256 element probes with a Vantage 256 system.  Probe
% %             % ID and other control signals (HVMux, temperature sense) will only
% %             % be active at connector 1.  Connector 2 is used only for element
% %             % signals.
% %             % First check for both connectors in same state
% %             if Presence(1) && ~Presence(2)
% %                 fprintf(2, 'VSX: Connection sensed at connector 1 but not connector 2.\n');
% %                 error('VSX: Both connectors must be in same state when using them together with Resource.Parameters.xconnector = 0.')
% %             elseif ~Presence(1) && Presence(2)
% %                 fprintf(2, 'VSX: Connection sensed at connector 2 but not connector 1.\n');
% %                 error('VSX: Both connectors must be in same state when using them together with Resource.Parameters.xconnector = 0.')
% %             elseif Presence(1) && Presence(2)
% %                 selectConnector = 0;
% %                 probeConnected = 1;
% %             elseif ~Presence(1) && ~Presence(2)
% %                 selectConnector = 0;
% %                 probeConnected = 0;
% %             end
% %         elseif isequal(Resource.SysConfig.UTAtype, [1 6 2 1])
% %             selectConnector = 0;
% %             probeConnected = Presence(2);
% %         end
%     else
%         % Resource.Parameters.Connector has length one, so select only that
%         % connector
    selectConnector = Resource.Parameters.Connector(1); % connector to be selected
    if selectConnector == 3
        % cannot sense connection at connector 3, so just assume connection
        % is present
        probeConnected = 1;
    else
        probeConnected = Presence(Resource.Parameters.Connector(1)); % current connection status
    end
%     end

    if probeConnected
        % a probe is connected so check if ID matches the value expected by
        % the script  (Note NDT UTA adapter has no ID support at any of the
        % three connectors, so if that adapter is present just act as if
        % the ID matches
        if isequal(Resource.SysConfig.UTAtype, [1 6 3 1]) || ID(max(selectConnector, 1)) == Trans.id
            % note that if we are selecting all connectors (selectConnector
            % = 0) then the ID will be taken from connector 1
            % Probe with the expected ID is connected, so select the
            % connector and proceed to run
            result = setSelectedConnector(Resource.Parameters.Connector);
            if ~strcmpi(result, 'SUCCESS')
                % report failure and exit
                error('VSX: Failed to select probe connector because "%s".\n', result);
            end
        else
            % a probe is connected with a different ID; check for 'custom'
            % transducer name to decide what to do next
            if strcmpi(Trans.name, 'custom')
                % Trans.name is 'custom' so user wants to ignore the actual
                % ID value read from the connected probe.  Before allowing
                % the script to run we have to see if the connected probe
                % uses HVMux chips; if so we will exit with an error
                % condition
                probeName = computeTrans(num2str(ID(selectConnector), '%04X'));
                if strcmp(probeName, 'Unknown')
                    hvMux = 0;
                else
                    % this is a recognized transducer so find HVMux status
                    hvMux = computeTrans(probeName, 'HVMux');
                end

                if hvMux
                    % an HVMux probe is connected so ignore the 'custom'
                    % script and exit with an error.
                    errorstr = ['VSX: Connected probe is ', probeName, ' using HVMux chips.  ID mismatch not allowed for HVMux probes.'];
                    error(errorstr);
                else
                    % not an hvMux probe, so go ahead and run the script
                    % normally, ignoring the ID mismatch
                    result = setSelectedConnector(Resource.Parameters.Connector);
                    if ~strcmpi(result, 'SUCCESS')
                        % report failure and exit
                        error('VSX: Failed to select probe connector because "%s".\n', result);
                    end
                    if Resource.Parameters.verbose>1
                        disp(['VSX status: ID of probe specified by script doesn''t match ID of connected probe (', probeName, '),']);
                        disp('but Trans.name in script is set to ''custom'' so ID mismatch is being ignored.');
                    end
                end
            else
                % Name is not 'custom' so just exit with an error
                % message:
                error('VSX: ID of probe specified by script doesn''t match ID of probe at connector %d.\n', selectConnector);
            end
        end
        % Select the connector for a valid ID.
    else
        % No probe detected at identified connector; determine if
        % fakeScanhead is set
        if isfield(Resource.Parameters,'fakeScanhead') && Resource.Parameters.fakeScanhead
            % fakeScanhead mode, so select connector and run
            if Resource.Parameters.verbose > 1
                % display status message
                disp('VSX Status: No probe has been detected at the specified connector, and ');
                disp('Resource.Parameters.fakeScanhead is set so proceeding to run the script');
                disp('on the HW system (this is live acquisition, not simulation mode).');
            end
            result = setSelectedConnector(Resource.Parameters.Connector);
            if ~strcmpi(result, 'SUCCESS')
                % report failure and exit
                error('VSX: Failed to select probe connector because "%s".\n', result);
            end
        else
            % nothing connected and fakeScanhead not set, so ask user if
            % they want to quit or automatically switch to fake mode.
            if Resource.Parameters.verbose
                disp(' ');
                disp('No probe detected at specified connector - ');
                disp('Do you want to run the script in fake scanhead mode for ');
                disp('live acquisition with nothing connected?');
                Reply = input('Enter Y or return to proceed, N to exit:', 's');
                if strcmpi(Reply, 'N')
                    % user says exit, so just return at this point
                    return
                end
                % user wants script to run, so set the flag to allow associated
                % logic to function normally
                Resource.Parameters.fakeScanhead = 1;
                result = setSelectedConnector(Resource.Parameters.Connector);
                if ~strcmpi(result, 'SUCCESS')
                    % report failure and exit
                    error('VSX: Failed to select probe connector because "%s".\n', result);
                end
            else
                % Resource.Parameters.verbose is false, so just exit with error message
                error('VSX: No probe detected at connector number %d.\n',selectConnector);
            end
        end
    end


    % If this scanhead is a HVMux scanhead, set the HVMux attributes.
    if isfield(Trans,'HVMux')
        if UTA.elBiasEna == 1 && Trans.elBias ~= 0
            error('VSX: HVMux-based probes cannot use Trans.elBias with this UTA module or SHI.');
        end
    end
    
    % Check for Element Bias use
    if Trans.elBias ~= 0
        % Element bias is requested; check status of UTA.elBiasEna
        switch UTA.elBiasEna
            case 0
                error('VSX: Trans.elBias can not be used with this UTA module.');
            case 1
                % check allowed range for HVMux power supply when used for
                % element bias
                if abs(Trans.elBias) < 10 || abs(Trans.elBias) > 100
                    error('VSX: Trans.elBias setting is outside the supported range of 10 to 100 Volts for HVMux power supply.');
                end
                Resource.VDAS.elBiasSel = 1;
            otherwise
                % note case 2 for UTA baseboard bias source will be added
                % in a future release
                error('VSX: Unrecognized UTA.elBiasEna value of "%d".\n', UTA.elBiasEna);
        end
    end
end

% Set up GUI; first check for user-defined replacement
if ~isfield(Resource.Parameters, 'GUI') || isempty(Resource.Parameters.GUI)
    Resource.Parameters.GUI = 'vsx_gui';
end
% now confirm the specified function actually exists and is on the path
if isempty(which(Resource.Parameters.GUI))
    error('VSX: The GUI function specified by ''Resource.Parameters.GUI'' could not be found.');
end
vsxGUIh = str2func(Resource.Parameters.GUI); % create the function call
vsxGUIh(); % and run it

% Set up DisplayWindow(s) for output image(s). DisplayWindow(s) can be created for displaying
% processed ImageBuffers or user data.  If no DisplayWindow structure is specified in the user script,
% this step is skipped.
if isfield(Resource,'DisplayWindow')
    for i = 1:size(Resource.DisplayWindow,2)
        % Set defaults for DisplayWindow structure
        if (~isfield(Resource.DisplayWindow(i),'Title'))||isempty(Resource.DisplayWindow(i).Title)
            Resource.DisplayWindow(i).Title = displayWindowTitle;
        end
        if (~isfield(Resource.DisplayWindow(i),'Position'))||isempty(Resource.DisplayWindow(i).Position)
            error('VSX: DisplayWindow(%d) has no ''Position'' attribute.',i);
        else
            imWidth = Resource.DisplayWindow(i).Position(3); % Position(3) is imWidth, not figure width.
            imHeight = Resource.DisplayWindow(i).Position(4); % Position(4) is imHeight, not figure height.
        end
        DisplayData = zeros(imHeight,imWidth,'uint8');
        if (~isfield(Resource.DisplayWindow(i),'Colormap'))||isempty(Resource.DisplayWindow(i).Colormap)
            Resource.DisplayWindow(1).Colormap = gray(256);  % default colormap is greyscale.
        end
        if (~isfield(Resource.DisplayWindow(i),'Type'))||isempty(Resource.DisplayWindow(i).Type)
            Resource.DisplayWindow(i).Type = 'Matlab';
        end
        switch Resource.DisplayWindow(i).Type
            case 'Matlab'  % DisplayWindow is Matlab figure window.
                % Create the figure window.
                Resource.DisplayWindow(i).figureHandle = figure( ...
                        'Name',Resource.DisplayWindow(i).Title,...
                        'NumberTitle','off',...
                        'Position',[Resource.DisplayWindow(i).Position(1), ... % left edge
                                    Resource.DisplayWindow(i).Position(2), ... % bottom
                                    imWidth + 100, imHeight + 150], ...            % width, height + border
                        'Colormap',Resource.DisplayWindow(i).Colormap, ...
                        'Visible','off');
                axes('Units','pixels','Position',[60,90,imWidth,imHeight]);
                set(gca, 'Units','normalized');  % restore normalized units for re-sizing window.
                if ~isfield(Resource.DisplayWindow(i),'splitPalette')||isempty(Resource.DisplayWindow(i).splitPalette)
                    Resource.DisplayWindow(i).splitPalette = 0;
                end
                if ~isfield(Resource.DisplayWindow(i),'pdelta')||isempty(Resource.DisplayWindow(i).pdelta)
                    Resource.DisplayWindow(i).pdelta = 0.5;
                end
                if ~isfield(Resource.DisplayWindow(i),'clrWindow')||isempty(Resource.DisplayWindow(i).clrWindow)
                    Resource.DisplayWindow(i).clrWindow = 0;
                end
                % Determine whether a 2D DisplayWindow is x,z (normal 2D scan), x,z or x,y (3D C-scan) oriented.
                % The DisplayWindow.ReferencePt is a point in the x,z plane for normal 2D scans, or a point in
                % the x,y,z volume for a 3D volume.  For 3D scans and DisplayWindow.mode = '2D', the
                % DisplayWindow.orientation attribute determines the orientation of the 2D slice.
                %  - check for ReferencePt specified and convert to 3 dimensions if necessary.
                if ~isfield(Resource.DisplayWindow(i),'ReferencePt')||isempty(Resource.DisplayWindow(i).ReferencePt)
                    Resource.DisplayWindow(i).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];
                end
                if length(Resource.DisplayWindow(i).ReferencePt)==2
                    Resource.DisplayWindow(i).ReferencePt(3) = Resource.DisplayWindow(i).ReferencePt(2);
                    Resource.DisplayWindow(i).ReferencePt(2) = 0;
                end
                if ~isfield(Resource.DisplayWindow(i),'mode')||isempty(Resource.DisplayWindow(i).mode)
                    Resource.DisplayWindow(i).mode = '2d';  % default mode to '2d'
                end
                if strcmp(Resource.DisplayWindow(i).mode,'2d')
                    if Trans.type < 2
                        % DisplayWindow is x,z oriented
                        vmin = Resource.DisplayWindow(i).ReferencePt(3); % vertical minimum is DisplayWindow.ReferencePt(3).
                        vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                        set(gca,'YDir','reverse');
                        xmin = Resource.DisplayWindow(i).ReferencePt(1);
                        xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
                   else
                        if ~isfield(Resource.DisplayWindow(i),'Orientation')
                            error(['VSX: Resource.DisplayWindow(%d).Orientation must be defined for a 2D Displaywindow',...
                                   'with Trans.type = 2\n'],i);
                        end
                        switch Resource.DisplayWindow(i).Orientation
                            case 'xz'
                                % DisplayWindow is x,z oriented
                                vmin = Resource.DisplayWindow(i).ReferencePt(3);
                                vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                                set(gca,'YDir','reverse');
                                xmin = Resource.DisplayWindow(i).ReferencePt(1);
                                xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
                            case 'yz'  % DisplayWindow is y,z oriented
                                vmin = Resource.DisplayWindow(i).ReferencePt(3);
                                vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                                set(gca,'YDir','reverse');
                                xmin = Resource.DisplayWindow(i).ReferencePt(2);
                                xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
                            case 'xy'  % DisplayWindow is x,y oriented
                                vmin = Resource.DisplayWindow(i).ReferencePt(2);
                                vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                                set(gca,'YDir','reverse');
                                xmin = Resource.DisplayWindow(i).ReferencePt(1);
                                xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
                        end
                    end
                else
                    error('VSX: Resource.DisplayWindow.modes other than ''2d'' not currently supported.');
                end
                % Set limits for x axis.
                scale = 1.0;  % default is no scaling (wavelengths)
                axisname = 'wavelengths';
                if isfield(Resource.DisplayWindow(i),'AxesUnits')&&~isempty(Resource.DisplayWindow(i).AxesUnits)
                    if strcmpi(Resource.DisplayWindow(i).AxesUnits,'mm')
                        axisname = 'mm';
                        if ~isfield(Resource.Parameters,'speedOfSound')||isempty(Resource.Parameters.speedOfSound)
                            scale = 1.54/Trans.frequency; % use default speed of sound.
                        else
                            scale = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
                        end
                    end
                end
                Resource.DisplayWindow(i).imageHandle = image('CData',DisplayData, ...
                        'XData',scale*[xmin,xmax], ...
                        'YData',scale*[vmin,vmax]);
                set(gca,'FontSize',12); xlabel(axisname,'FontSize',14);
                axis equal tight;
                drawnow
                set(Resource.DisplayWindow(i).figureHandle,'Visible','on');
            case 'Verasonics'  % DisplayWindow is Verasonics Java window.

                % Import and initialize the Java libaries that are needed for the Verasonics viewer.
                import com.verasonics.viewer.ui.components.*
                import com.verasonics.viewer.ui.panels.*
                import com.verasonics.viewer.ui.*
                import com.verasonics.viewer.*
                import com.verasonics.vantage.image.*
                import com.verasonics.vantage.image.events.*

                % Uncomment the line below to print image viewer events (enable debug mode).
                % MatlabVantageWindow.setDebugMode(true);

                % Define the background color of the viewer window.
                viewerBgColor = '0xCDCDCD';

                % Create the Verasonics Window.
                vantageWindow = MatlabVantageWindow.create(Attr('xPos',Resource.DisplayWindow(i).Position(1)), ... % left edge
                                                           Attr('yPos',Resource.DisplayWindow(i).Position(2)), ... % bottom
                                                           Attr('width',imWidth + VantageWindow.PAD_WIDTH), ...    % width of image + padding
                                                           Attr('height',imHeight + VantageWindow.PAD_HEIGHT), ...  % height of image + padding
                                                           Attr('resizeable', true), ...
                                                           Attr('bkgColor', viewerBgColor), ...
                                                           Attr('title',Resource.DisplayWindow(i).Title));

                % Get the Unique ID of this window. Used to communicate with this specific window.
                Resource.DisplayWindow(i).figureHandle = vantageWindow;

                % Create the image viewer that is displayed in the window.
                imageViewer = ImageViewerPanel(Attr('label',Resource.DisplayWindow(i).Title), ...
                                               Attr('showAxisAnnotation', true), ...
                                               Attr('bkgColor', viewerBgColor), ...
                                               Attr('scaleToWindowSize', true), ...
                                               Attr('acqRateMethod', 'synchConstantRate30Hz'));

                 % Add the image viewer to the window.
                vantageWindow.addComponent(imageViewer, ... % component
				                           0, 0, ... % X|Y grid placement
                                           Attr('anchor', 'north'), ...
                                           Attr('weightX', 1.0), ...
				                           Attr('weightY', 1.0));

				% Add the menubar to the image viewer window.
                imageViewerMenu = MatlabImageViewerMenu(vantageWindow, imageViewer);
                vantageWindow.addMenuBar(Attr('menuBarObject', imageViewerMenu));

                % Get the Unique ID of this image viewer. Used to send events to this specific image viewer.
                Resource.DisplayWindow(i).imageHandle = imageViewer.uid;

                % Set axes limits base on Resource.DisplayWindow(i).Orientation
                if Trans.type < 2
                    % DisplayWindow is x,z oriented
                    vmin = Resource.DisplayWindow(i).ReferencePt(3); % vertical minimum is DisplayWindow.ReferencePt(3).
                    vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                    xmin = Resource.DisplayWindow(i).ReferencePt(1);
                    xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
               else
                    if ~isfield(Resource.DisplayWindow(i),'Orientation')
                        error(['VSX: Resource.DisplayWindow(%d).Orientation must be defined for a 2D Displaywindow',...
                               'with Trans.type = 2\n'],i);
                    end
                    switch Resource.DisplayWindow(i).Orientation
                        case 'xz'
                            % DisplayWindow is x,z oriented
                            vmin = Resource.DisplayWindow(i).ReferencePt(3);
                            vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                            xmin = Resource.DisplayWindow(i).ReferencePt(1);
                            xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
                        case 'yz'  % DisplayWindow is y,z oriented
                            vmin = Resource.DisplayWindow(i).ReferencePt(3);
                            vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                            xmin = Resource.DisplayWindow(i).ReferencePt(2);
                            xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
                        case 'xy'  % DisplayWindow is x,y oriented
                            vmin = Resource.DisplayWindow(i).ReferencePt(2);
                            vmax = vmin + imHeight*Resource.DisplayWindow(i).pdelta;
                            xmin = Resource.DisplayWindow(i).ReferencePt(1);
                            xmax = xmin + imWidth*Resource.DisplayWindow(i).pdelta;
                    end
                end
                % Set limits according to AxesUnits.
                scale = 1.0;  % default is no scaling (wavelengths)
                axisname = 'wavelengths';
                unitsType = SpatialUnitsType.wavelengths;
                if isfield(Resource.DisplayWindow(i),'AxesUnits')&&~isempty(Resource.DisplayWindow(i).AxesUnits)
                    if strcmpi(Resource.DisplayWindow(i).AxesUnits,'mm')
                        axisname = 'mm';
                        unitsType = SpatialUnitsType.millimeters;
                        if ~isfield(Resource.Parameters,'speedOfSound')||isempty(Resource.Parameters.speedOfSound)
                            scale = 1.54/Trans.frequency; % use default speed of sound.
                        else
                            scale = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
                        end
                    end
                end

                % Generate an event to update the colormap.
                colorMapEvent = ColorMapEvent(256, ColorMapType.planarNoramlizedDoubles, Resource.DisplayWindow(i).Colormap);
                VantageImageEvent.generateColorMapEvent(imageViewer.uid, colorMapEvent);

                % Generate an event that the spatial units have changed.
                spatialUnitsEvent = SpatialUnitsEvent(unitsType, scale);
                VantageImageEvent.generateSpatialUnitsEvent(imageViewer.uid, spatialUnitsEvent);

                % Generate an event that the spatial position has changed.
                spatialPositionEvent = SpatialPositionEvent(xmin, xmax, vmin, vmax);
                VantageImageEvent.generateSpatialPositionEvent(imageViewer.uid, spatialPositionEvent);

            otherwise
                error('VSX: Unrecognized Type for DisplayWindow(%d).',i);
        end
    end

    clear imWidth imHeight
    % Capture Resource.DisplayWindow(1).pdelta for use by zoom functions.
    orgDispPdelta = Resource.DisplayWindow(1).pdelta;

end

% ***** Initialize GUI and add UI objects. *****
% - Set state of simulate and rcvdataloop buttons (initial state is 0).
set(findobj('String','Rcv Data Loop'),'Value',rloopButton);
set(findobj('String','Simulate'),'Value',simButton);

% Make UI window the current figure, unless the caller has requested to hide it.
if((true == exist('Mcr_GuiHide', 'var')) && (1 == Mcr_GuiHide))
  % Caller has requested that we do NOT show the GUI window.
  if Resource.Parameters.verbose>1
      display('NOTE:  Caller has requested that the VSX GUI window be hidden.  Use ctrl-c to abort VSX if necessary.');
  end
else
  % OK, DO make the GUI window active.
  set(0, 'CurrentFigure', findobj('tag','UI'));
end

% - If .mat file contains a UI object, add it to UI window.
if exist('UI','var')
    if strcmp('vsx_gui',Resource.Parameters.GUI)
        % Get GUI background color.
        f = findobj('tag','UI');
        bkgrnd = get(f,'Color');
        % Define available user control labels and their positions.
        UserID = {'UserA1','UserA2','UserB1','UserB2','UserB3','UserB4','UserB5','UserB6','UserB7','UserB8',...
                  'UserC1','UserC2','UserC3','UserC4','UserC5','UserC6','UserC7','UserC8'};
        UserPos = [UIPos(2,:,1); UIPos(3,:,1);...
                   UIPos(2,:,2); UIPos(3,:,2); UIPos(4,:,2); UIPos(5,:,2); UIPos(6,:,2); UIPos(7,:,2); UIPos(8,:,2); UIPos(9,:,2);...
                   UIPos(2,:,3); UIPos(3,:,3); UIPos(4,:,3); UIPos(5,:,3); UIPos(6,:,3); UIPos(7,:,3); UIPos(8,:,3); UIPos(9,:,3)];
    end
    for i = 1:size(UI,2)
        if (isfield(UI(i),'Statement')&&~isempty(UI(i).Statement))
            eval(UI(i).Statement);
        end
        if isfield(UI(i),'Control') && ~isempty(UI(i).Control)
            if strcmp('vsx_gui',Resource.Parameters.GUI)
                if strncmp(UI(i).Control{1},'User',4) % Do we have a UserXX control?
                    L = strcmp('Style',UI(i).Control);
                    j = find(L,1);
                    if isempty(j), error('VSX: No ''Style'' defined for user UI(%d).Control.\n',i); end
                    UIStyle = UI(i).Control{j+1};
                    L = strcmp(UI(i).Control{1},UserID);
                    j = find(L,1);
                    if isempty(j), error('VSX: No definition exists for UI(%d).Control at %s.\n',i,UI(i).Control{1}); end
                    switch UIStyle
                        case 'VsSlider'
                            L = strcmpi('Label',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), Txt = UI(i).Control{k+1}; else Txt = UserID(j); end
                            L = strcmpi('SliderMinMaxVal',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), MMV = UI(i).Control{k+1}; else MMV = [1,100,1]; end
                            L = strcmpi('SliderStep',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), SS = UI(i).Control{k+1}; else SS = [0.01,0.1]; end
                            L = strcmpi('ValueFormat',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), VF = UI(i).Control{k+1}; else VF = '%3.0f'; end
                            UserT = {f,'Style','text','String',Txt,'Units','normalized',...
                                'Position',[(UserPos(j,:)+SG.TO) SG.TS],...
                                'FontUnits','normalized','FontSize',SG.TF,'FontWeight','bold'};
                            UserS = {f,'Style','slider','Min',MMV(1),'Max',MMV(2),'Value',MMV(3),...
                                'SliderStep',SS,'Units','normalized','Position',[UserPos(j,:)+SG.SO,SG.SS],...
                                'BackgroundColor',bkgrnd-0.05,'Tag',[UserID{j} 'Slider'],...
                                'Callback',{str2func([UserID{j} 'Callback'])}};
                            UserE = {f,'Style','edit','String',num2str(MMV(3),VF),'UserData',VF,...
                                'Units','normalized','Position',[(UserPos(j,:)+SG.EO) SG.ES],...
                                'BackgroundColor',bkgrnd+0.1,'Tag',[UserID{j} 'Edit'],...
                                'Callback',{str2func([UserID{j} 'Callback'])}};
                            UI(i).handle = [0,0,0];
                            UI(i).handle(1) = uicontrol(UserT{:});
                            UI(i).handle(2) = uicontrol(UserS{:});
                            UI(i).handle(3) = uicontrol(UserE{:});
                            clear VF
                            % Define UserXX.Callback preamble.
                            UserCBPre = {...
                                [UserID{j} 'Callback(hObject, eventdata);'],...
                                ' ',...
                                'Cntrl = get(hObject, ''Style'');',...
                                'if strcmp(Cntrl,''slider'')',...
                                '    UIValue = get(hObject,''Value'');',...
                                ['    h = findobj(''Tag'', ''' [UserID{j} 'Edit'] ''');'],...
                                '    set(h,''String'',num2str(UIValue,get(h,''UserData'')));',...
                                'else',...
                                '    UIValue = str2num(get(hObject,''String''));',...
                                ['    h = findobj(''Tag'', ''' [UserID{j} 'Slider'] ''');'],...
                                '    mx = get(h,''Max'');',...
                                '    mn = get(h,''Min'');',...
                                '    if (UIValue > mx)',...
                                '        UIValue = mx;',...
                                '        set(hObject,''String'',num2str(UIValue,get(hObject,''UserData'')));',...
                                '    end',...
                                '    if (UIValue < mn)',...
                                '        UIValue = mn;',...
                                '        set(hObject,''String'',num2str(UIValue,get(hObject,''UserData'')));',...
                                '    end',...
                                '    set(h,''Value'',UIValue);',...
                                'end'};
                        case 'VsPushButton'
                            L = strcmpi('Label',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), Txt = UI(i).Control{k+1}; else Txt = UserID(j); end
                            UI(i).handle = uicontrol(f,'Style','pushbutton','String',Txt,....
                                'Units','normalized','Position',[UserPos(j,:)+PB.BO,PB.BS],...
                                'FontUnits','normalized','FontSize',PB.FS,...
                                'BackgroundColor',bkgrnd+0.05,...
                                'Callback',{str2func([UserID{j} 'Callback'])});
                            % Define UserXX.Callback preamble.
                            UserCBPre = {...
                                [UserID{j} 'Callback(hObject, eventdata);'],...
                                ' '};
                        case 'VsToggleButton'
                            L = strcmpi('Label',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), Txt = UI(i).Control{k+1}; else Txt = UserID(j); end
                            UI(i).handle = uicontrol(f,'Style','togglebutton','String',Txt,...
                                'Units','normalized','Position',[UserPos(j,:)+PB.BO,PB.BS],...
                                'FontUnits','normalized',...    %'FontSize',PB.FS,...
                                'BackgroundColor',bkgrnd+0.05,'tag',[UserID{j} 'ToggleButton'],...
                                'Callback',{str2func([UserID{j} 'Callback'])});
                            % Define UserXX.Callback preamble.
                            UserCBPre = {...
                                [UserID{j} 'Callback(hObject, eventdata);'],...
                                ' ',...
                                'button_state = get(hObject,''Value'');',...
                                'if button_state == get(hObject,''Max'')',...
                                '   % Toggle button is pressed.',...
                                '   UIState = 1;',...
                                'elseif button_state == get(hObject,''Min'')',...
                                '   % Toggle button is not pressed.',...
                                '   UIState = 0;',...
                                'end'};
                        case 'VsButtonGroup'
                            L = strcmpi('Title',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), Txt = UI(i).Control{k+1}; else Txt = UserID(j); end
                            L = strcmpi('NumButtons',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), nbuttons = UI(i).Control{k+1}; else nbuttons = 2; end
                            L = strcmpi('Labels',UI(i).Control);
                            k = find(L,1);
                            if ~isempty(k), Labels = UI(i).Control{k+1};
                            else for k=1:nbuttons, Labels(k) = num2str(k); end
                            end
                            % Modify sizes for no. of buttons.
                            BG.BGS(2) = (nbuttons+1)*BG.BI;
                            BG.BGO(2) = 0.12 - BG.BGS(2);
                            BG.BS(2) = 1.08/(nbuttons+1);
                            BG.BO(2) = 0.08/(nbuttons-1);%1-2*BG.BS(2);
                            BG.TFS = 1/(nbuttons+2.5);
                            UI(i).handle = zeros(1,nbuttons);
                            UI(i).handle(1) = uibuttongroup('Title',Txt,'Units','normalized',...
                                'Position',[UserPos(j,:)+BG.BGO,BG.BGS],...
                                'FontUnits','normalized','FontSize',BG.TFS,'FontWeight','bold',...
                                'BackgroundColor',bkgrnd,'SelectionChangeFcn',{str2func([UserID{j} 'Callback'])});
                            for k=1:nbuttons
                                UI(i).handle(k+1) = uicontrol('Style','radiobutton','Parent',UI(i).handle(1),...
                                    'Units','normalized',...
                                    'Position',[BG.BO(1),BG.BO(2)+(nbuttons-k)*(1/nbuttons),BG.BS],...
                                    'FontUnits','normalized','FontSize',BG.BFS,'String',Labels(k),...
                                    'Tag',[UserID{j} 'RadioButton' num2str(k)]);
                            end
                            % Define UserXX.Callback preamble.
                            UserCBPre = {...
                                [UserID{j} 'Callback(hObject, eventdata);'],...
                                ' ',...
                                'S = get(eventdata.NewValue,''Tag'');',...
                                'UIState = str2double(S(18));',...
                                };
                        otherwise
                            error('VSX: Unrecognized style for user UI(%d).Control.',i);
                    end
                else
                    UI(i).handle = uicontrol(UI(i).Control{:});
                end
            else
                if strcmp('vsx_guiPC',Resource.Parameters.GUI)
                % Adjust Windows UIControl parameters to match Unix values
                    L = strcmpi('Style',UI(i).Control);
                    index = find(L,1);
                    if ~isempty(index)
                        switch UI(i).Control{index+1}
                            case 'text'
                                L = strcmpi('Position',UI(i).Control);
                                index = find(L,1);
                                if ~isempty(index)
                                    UI(i).Control{index+1}(2) = UI(i).Control{index+1}(2)-4;
                                end
                            case 'slider'
                                L = strcmpi('Position',UI(i).Control);
                                index = find(L,1);
                                if ~isempty(index)
                                    UI(i).Control{index+1}(4) = UI(i).Control{index+1}(4)-10;
                                end
                            case 'edit'
                                L = strcmpi('Position',UI(i).Control);
                                index = find(L,1);
                                if ~isempty(index)
                                    UI(i).Control{index+1}(1) = UI(i).Control{index+1}(1)+5;
                                    UI(i).Control{index+1}(2) = UI(i).Control{index+1}(2)-12;
                                end
                        end
                    end
                end
                UI(i).handle = uicontrol(UI(i).Control{:});
            end
        end
        if (isfield(UI(i),'Callback')&&~isempty(UI(i).Callback))
            % Decode callback function.
            %
            % WARNING:  MCR cannot use dynamically generated .m scripts.  They
            % must exist at the time MATLAB Compiler was run.
            %
            % Thus, use "ismcc"/"isdeployed" to avoid dynamically generated .m
            % files.  For Verasonics SetUpXxx.m calls, this means the
            % Verasonics UI.Callback entries that use tempdir dynamically
            % generated scripts.
            if(false == isdeployed)
                % OK, we are NOT under MCR, so DO make these calls.
                if (isfield(UI(i),'Control')&&~isempty(UI(i).Control))&&strncmp(UI(i).Control{1},'User',4) % Do we have a UserXX control?
                    fname = textscan(UserCBPre{1},'%s %*[^\n]', 'Delimiter','(');
                    fid = fopen([tempdir,fname{1}{1},'.m'], 'w');
                    fprintf(fid,'function %s\n', UserCBPre{1});
                    for j = 2:size(UserCBPre,2)
                        fprintf(fid,'%s\n', UserCBPre{j});
                    end
                    for j = 1:size(UI(i).Callback,2)
                        fprintf(fid,'%s\n', UI(i).Callback{j});
                    end
                else
                    % Check for older format that required callback name on first line.
                    S = deblank(UI(i).Callback{1});
                    S = S((length(S)-1):length(S));  % if older format, S will get '.m'
                    if strcmp(S,'.m') % old format
                        fid = fopen([tempdir,UI(i).Callback{1}], 'w');
                        for j = 2:size(UI(i).Callback,2)
                            fprintf(fid,'%s\n', UI(i).Callback{j});
                        end
                    else % new format that has function prototype as first line.
                        n = strfind(UI(i).Callback{1},'=');
                        if ~isempty(n)
                            fname = textscan(UI(i).Callback{1},'%*s %s %*[^\n]','Delimiter',{'=','('});
                        else
                            fname = textscan(UI(i).Callback{1},'%s %*[^\n]', 'Delimiter','(');
                        end
                        fid = fopen([tempdir,fname{1}{1},'.m'], 'w');
                        fprintf(fid,'function %s\n', UI(i).Callback{1});
                        for j = 2:size(UI(i).Callback,2)
                            fprintf(fid,'%s\n', UI(i).Callback{j});
                        end
                    end
                end
                status = fclose(fid);
                addpath(tempdir);
                rehash path
            end
        end
    end
end


%% ***** HIFU Limits Check *****
% this check is now made in the initialize function of runAcq, after
% VsUpdate(TW) and VsUpdate(TX) have been completed.

%% ***** Start of loop for continuously calling 'runAcq'. *****

% TEST: return here to examine workspace without calling runAcq
% return

% 'runAcq' is called with the single input, 'Control', which is a two attribute structure
%    array where the first attribute, 'Command', is a command string,  and the second attribute,
%    'Parameters', is a cell array whose interpretation depends on the command given.
% Valid commands are:
%    Control.Command =
%      'set' - set an attribute of an existing object (structure) and return.
%          Control.Parameters = {'Object', objectNumber, 'attribute', value, 'attribute', value, ...}
%      'set&Run' - set an attribute of an existing object (structure) and run sequence.
%          Control.Parameters = {'Object', objectNumber, 'attribute', value, 'attribute', value, ...}
%      'update&Run' - update the entire object(s) by re-reading from Matlab, then run sequence.
%          Control.Parameters = {'InterBuffer','ImageBuffer','DisplayWindow',
%                                'Parameters', 'Trans','HVMux','Media','PData',
%                                'TW','TX','Receive','ReceiveProfiles','TGC',
%                                'Recon','Process','SeqControl','Event'}*
%                                * provide one or more of this list of objects
%      'copyBuffers' - Copies back to Matlab IQData,ImgData,ImgDataP and LogData
%          buffers, without running the Event sequence.
%      'setBFlag' - set the BFlag used for triggering the conditional branch sequence control.
%      'imageDisplay' - display an image using the specified Image parameters.
%          Control.Parameters = {'Image attribute, value, ...}
%      'debug' = turns on or off debug output.
if ~exist('Control', 'var')
    Control.Command = [];  % Set Control.Command to 'empty' for no control action.
end
%Control.Command = 'debug';  % Control.Command to turn on or off debug output.
%Control.Parameters = {'on'};
action = 'none';  % First main loop action must be 'none'.
exit = 0;   % exit will get set to 1 when UI window is closed.
freeze = 0; % freeze will get set to 1 when 'Freeze' togglebutton is pressed.
sequencePeriod = 1;   % set initial frame period to 1 second
initialized = 0;  % will be set to one by runAcq after hardware initialization
tElapsed = 0;

%% final initialization step: calls to VsUpdate and TXEventCheck
% note these initialization calls must be made in the order listed here
% an error from any of these function calls will result in VSX exiting
% with an error
updateh('TW');
updateh('TX');
updateh('Receive');
updateh('SeqControl');
if VDASupdates; TXEventCheckh(); end


% Check to see if user only wants runAcq to initialize.
if Resource.Parameters.initializeOnly > 0
    runAcq(Control);
    return
end

%% Run sequence.
while exit==0
    drawnow
    if freeze == 1
        % If running with the hardware, stop the hardware sequencer and wait for freeze
        % togglebutton to be pressed off.
        if (VDAS==1), Result = sequenceStop; end
        % We are about to enter a wait-for-unfreeze.  During that time, no commands besides
        % callbacks are available.  For callers that require the ability to send a command during
        % the time VSX is frozen, we implement a MATLAB Timer.
        %
        % If the caller sets the variable Mcr_FreezeTimerFunction to the name of their callback
        % function, then a timer will be created to periodically call that callback, thus allowing
        % the injection of commands while frozen, such as unfreeze.  (This would be the case for a
        % caller who has hidden the VSX GUI control panel, and thus the unfreeze button, with their
        % own custom GUI control panel that uses only "pull" technology rather than VSX GUI control
        % panel's "push" technology.  "Pull" here means that case where VSX is running separately
        % from a custom control panel and VSX is asking the custom control panel for any commands
        % to execute.
        %
        % NOTE:  Because the timer works within the MATLAB single-threaded
        % environment, it cannot guarantee execution times or execution rates.
        if((true == exist('Mcr_FreezeTimerFunction', 'var')) && (0 == strcmp('', Mcr_FreezeTimerFunction)))
          % OK, we DO have a callback hook to call from a Timer.
          %
          % We want to call the user's callback SCRIPT, but Timers want to
          % call a FUNCTION that accepts (object,event).  So, we make an
          % anonymous function here that calls the user's script.
          functionHandle = @(object, event)evalin('base', Mcr_FreezeTimerFunction);
          freezeTimer = timer ...
          ( ...
            'BusyMode',      'queue', ...        % So communication is not lost
            'ExecutionMode', 'fixedRate', ...
            'StartDelay',    0, ...              % Call callback immediately
            'Period',        0.25, ...           % Call again every 1/4 seconds
            'TimerFcn',      functionHandle ... % User's callback
          );

          start(freezeTimer);
        end

        % Wait for unfreeze
        waitfor(findobj('String','Freeze'),'Value',0);
        if exit==1 % exit might get set while in freeze.
            break;
        end

        % If we used a timer above, stop it and delete it.
        if(true == exist('freezeTimer', 'var'))
          % OK, we DO have a timer to stop and delete.
          stop(freezeTimer);
          delete(freezeTimer);
          clear freezeTimer;
        end

        % Reset startEvent, in case sequence contains return-to-Matlab points prior to the
        % end of sequence, where 'freeze' may have been pressed.  Sequence will be re-started
        % after updating startEvent.
        if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
        Control(n).Command = 'set&Run';
        Control(n).Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
    end
    % The following switch statement implements actions set by the default GUI controls.
    switch action
        case 'none'
        case 'tgc'
            nTGC = 1;
            if length(TGC) > 1
                TGCnum = findobj('Tag','TGCnum');
                nTGC = get(TGCnum,'Value');
            end
            TGC(nTGC).CntrlPts = [tgc1,tgc2,tgc3,tgc4,tgc5,tgc6,tgc7,tgc8];
            TGC(nTGC).Waveform = computeTGCWaveform(TGC(nTGC));
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'set&Run';
            Control(n).Parameters = {'TGC', nTGC, 'Waveform', TGC(nTGC).Waveform};
        case 'zoomin'
            % 'zoomin' keeps the size of the DisplayWindow constant, but decreases both
            % DisplayWindow.pdelta and PData(1).pdelta. The DisplayWindow.refPt also is changed to
            % zoom from the center. Note: 'zoomin' only works on DisplayWindow 1, which is assumed
            % to be the main image display window.
            cx = Resource.DisplayWindow(1).ReferencePt(1) + Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(3)-1)/2;
            cz = Resource.DisplayWindow(1).ReferencePt(3) + Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(4)-1)/2;
            Resource.DisplayWindow(1).pdelta = Resource.DisplayWindow(1).pdelta * 0.8;
            if Resource.DisplayWindow(1).pdelta < orgDispPdelta * 0.2621
                Resource.DisplayWindow(1).pdelta = orgDispPdelta * 0.2621;
            end
            Resource.DisplayWindow(1).ReferencePt(1) = cx - Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(3)-1)/2;
            Resource.DisplayWindow(1).ReferencePt(3) = cz - Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(4)-1)/2;
            PData(1).Origin(1) = Resource.DisplayWindow(1).ReferencePt(1);
            PData(1).Origin(3) = Resource.DisplayWindow(1).ReferencePt(3);
            PData(1).PDelta(1) = PData(1).PDelta(1) * 0.8;
            if PData(1).PDelta(1) < orgPdeltaX * 0.2621, PData(1).PDelta(1) = orgPdeltaX * 0.2621; end
            PData(1).PDelta(3) = PData(1).PDelta(3) * 0.8;
            if PData(1).PDelta(3) < orgPdeltaZ * 0.2621, PData(1).PDelta(3) = orgPdeltaZ * 0.2621; end
            axesHandle = get(Resource.DisplayWindow(1).figureHandle,'CurrentAxes');
            cla(axesHandle); % this clears any objects attached to the old axes.
            ymin = Resource.DisplayWindow(1).ReferencePt(3);
            ymax = ymin + Resource.DisplayWindow(1).Position(4)*Resource.DisplayWindow(1).pdelta;
            xmin = Resource.DisplayWindow(1).ReferencePt(1);
            xmax = xmin + Resource.DisplayWindow(1).Position(3)*Resource.DisplayWindow(1).pdelta;
            DisplayData = zeros(Resource.DisplayWindow(1).Position(4),Resource.DisplayWindow(1).Position(3),'uint8');
            Resource.DisplayWindow(1).imageHandle = image(DisplayData,...
                    'Parent',axesHandle, ...
                    'Xdata',scale*[xmin,xmax], ...
                    'YData',scale*[ymin,ymax]);
            Resource.DisplayWindow(1).clrWindow = 1;
            Resource.ImageBuffer(1).lastFrame = 0; % this clears the ImageBuffer
            PData(1).Region = computeRegions(PData(1));
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'update&Run';
            Control(n).Parameters = {'PData','Recon','DisplayWindow','ImageBuffer'};
        case 'zoomout'
            cx = Resource.DisplayWindow(1).ReferencePt(1) + Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(3)-1)/2;
            cz = Resource.DisplayWindow(1).ReferencePt(3) + Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(4)-1)/2;
            Resource.DisplayWindow(1).pdelta = Resource.DisplayWindow(1).pdelta * 1.25;
            if Resource.DisplayWindow(1).pdelta > orgDispPdelta * 3.8147
                Resource.DisplayWindow(1).pdelta = orgDispPdelta * 3.8147;
            end
            Resource.DisplayWindow(1).ReferencePt(1) = cx - Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(3)-1)/2;
            Resource.DisplayWindow(1).ReferencePt(3) = cz - Resource.DisplayWindow(1).pdelta*...
                (Resource.DisplayWindow(1).Position(4)-1)/2;
            PData(1).Origin(1) = Resource.DisplayWindow(1).ReferencePt(1);
            PData(1).Origin(3) = Resource.DisplayWindow(1).ReferencePt(3);
            PData(1).PDelta(1) = PData(1).PDelta(1) * 1.25;
            if PData(1).PDelta(1) > orgPdeltaX * 3.8147, PData(1).PDelta(1) = orgPdeltaX * 3.8147; end
            PData(1).PDelta(3) = PData(1).PDelta(3) * 1.25;
            if PData(1).PDelta(3) > orgPdeltaZ * 3.8147, PData(1).PDelta(3) = orgPdeltaZ * 3.8147; end
            axesHandle = get(Resource.DisplayWindow(1).figureHandle,'CurrentAxes');
            cla(axesHandle); % this clears any objects attached to the old axes.
            ymin = Resource.DisplayWindow(1).ReferencePt(3);
            ymax = ymin + Resource.DisplayWindow(1).Position(4)*Resource.DisplayWindow(1).pdelta;
            xmin = Resource.DisplayWindow(1).ReferencePt(1);
            xmax = xmin + Resource.DisplayWindow(1).Position(3)*Resource.DisplayWindow(1).pdelta;
            DisplayData = zeros(Resource.DisplayWindow(1).Position(4),Resource.DisplayWindow(1).Position(3),'uint8');
            Resource.DisplayWindow(1).imageHandle = image(DisplayData,...
                    'Parent',axesHandle, ...
                    'Xdata',scale*[xmin,xmax], ...
                    'YData',scale*[ymin,ymax]);
            Resource.DisplayWindow(1).clrWindow = 1;
            Resource.ImageBuffer(1).lastFrame = 0; % this clears the ImageBuffer
            PData(1).Region = computeRegions(PData(1));
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'update&Run';
            Control(n).Parameters = {'PData','Recon','DisplayWindow','ImageBuffer'};
        case 'panlft'
            % pan 10% of image width.
            PData(1).Origin(1) = PData(1).Origin(1) - 0.1*PData(1).PDelta(1)*PData(1).Size(2);
            Resource.DisplayWindow(1).ReferencePt(1) = PData(1).Origin(1);
            Resource.DisplayWindow(1).clrWindow = 1;
            Resource.ImageBuffer(1).lastFrame = 0;
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'update&Run';
            Control(n).Parameters = {'DisplayWindow'};
        case 'panrt'
            % pan 10% of image width.
            PData(1).Origin(1) = PData(1).Origin(1) + 0.1*PData(1).PDelta(1)*PData(1).Size(2);
            Resource.DisplayWindow(1).ReferencePt(1) = PData(1).Origin(1);
            Resource.DisplayWindow(1).clrWindow = 1;
            Resource.ImageBuffer(1).lastFrame = 0;
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'update&Run';
            Control(n).Parameters = {'DisplayWindow'};
        case 'panup'
            % pan 10% of image height.
            PData(1).Origin(3) = PData(1).Origin(3) - 0.1*PData(1).PDelta(3)*PData(1).Size(1);
            Resource.DisplayWindow(1).ReferencePt(3) = PData(1).Origin(3);
            Resource.DisplayWindow(1).clrWindow = 1;
            Resource.ImageBuffer(1).lastFrame = 0;
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'update&Run';
            Control(n).Parameters = {'DisplayWindow'};
        case 'pandn'
            % pan 10% of image height.
            PData(1).Origin(3) = PData(1).Origin(3) + 0.1*PData(1).PDelta(3)*PData(1).Size(1);
            Resource.DisplayWindow(1).ReferencePt(3) = PData(1).Origin(3);
            Resource.DisplayWindow(1).clrWindow = 1;
            Resource.ImageBuffer(1).lastFrame = 0;
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'update&Run';
            Control(n).Parameters = {'DisplayWindow'};
        case 'persist'
            % Set the value of persistence for the first Process structure with class = 'image, method =
            %  'imageDisplay' and 'displayWindow' = 1;
            flag = 0;
            if exist('Process','var')
                % Find the first 'Image/imageDisplay' Process structure for displayWindow 1
                for i = 1:size(Process,2)
                    if (strcmp(Process(i).classname,'Image'))&&(strcmp(Process(i).method,'imageDisplay'))
                        for j = 1:2:size(Process(i).Parameters,2)
                            if (strcmp(Process(i).Parameters{j},'displayWindow'))&&(Process(i).Parameters{j+1}==1)
                                for k = 1:2:size(Process(i).Parameters,2)
                                    if strcmp(Process(i).Parameters{k},'persistLevel')
                                        Process(i).Parameters{k+1} = persist;
                                        if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                                        Control(n).Command = 'set&Run';
                                        Control(n).Parameters = {'Process',i,...
                                                              'persistLevel',persist};
                                        flag = 1;
                                    end
                                end
                            end
                        end
                    end
                    if flag==1, break, end
                end
            end
        case 'pgain'
            % change processing gain factor for 1st 'image' Process structure with method "imageDisplay".
            for i = 1:size(Process,2);
                if (strcmp(Process(i).classname,'Image'))&&(strcmp(Process(i).method,'imageDisplay'))
                    break;
                end
            end
            if i <= size(Process,2)
                if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                Control(n).Command = 'set&Run';
                Control(n).Parameters = {'Process',i,...
                                         'pgain',pgain};
            end
        case 'speed'
            % change speed of sound correction factor.
            Resource.Parameters.speedCorrectionFactor = speedCorrect;
            if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
            Control(n).Command = 'set&Run';
            Control(n).Parameters = {'Parameters',1,...
                                  'speedCorrectionFactor',speedCorrect};
        case 'displayChange'
            for i = 1:numel(Resource.DisplayWindow)
                % adjust Matlab display window for size, position or aspect ratio change.
                if ~isempty(Resource.DisplayWindow(i).Type)&&strcmp(Resource.DisplayWindow(i).Type,'Matlab')
                    if i > 1, continue; end  % Don't modify any Matlab windows other than first.
                    axesHandle = get(Resource.DisplayWindow(i).figureHandle,'CurrentAxes');
                    cla(axesHandle); % this clears any objects attached to the old axes.
                    oldPos = get(Resource.DisplayWindow(i).figureHandle,'Position');
                    set(Resource.DisplayWindow(i).figureHandle, ...
                        'Position',[oldPos(1), ... % left edge
                                    oldPos(2), ... % bottom
                                    Resource.DisplayWindow(i).Position(3) + 100, ... % width
                                    Resource.DisplayWindow(i).Position(4) + 150]);    % height
                    set(axesHandle, 'Units','pixels', ...
                                    'Position',[60,90,Resource.DisplayWindow(i).Position(3),Resource.DisplayWindow(i).Position(4)]);
                    set(axesHandle, 'Units','normalized');  % restore normalized units for re-sizing window.
                    ymin = Resource.DisplayWindow(i).ReferencePt(3);
                    ymax = ymin + Resource.DisplayWindow(i).Position(4)*Resource.DisplayWindow(i).pdelta;
                    xmin = Resource.DisplayWindow(i).ReferencePt(1);
                    xmax = xmin + Resource.DisplayWindow(i).Position(3)*Resource.DisplayWindow(i).pdelta;
                    DisplayData = zeros(Resource.DisplayWindow(i).Position(4),Resource.DisplayWindow(i).Position(3),'uint8');
                    Resource.DisplayWindow(i).imageHandle = image(DisplayData,...
                            'Parent',axesHandle, ...
                            'Xdata',scale*[xmin,xmax], ...
                            'YData',scale*[ymin,ymax]);
                    Resource.DisplayWindow(i).clrWindow = 1;
                    Resource.ImageBuffer(i).lastFrame = 0; % this clears the ImageBuffer
                    axis(axesHandle, 'equal', 'tight');
                    set(axesHandle,'FontSize',12); xlabel(axesHandle,axisname,'FontSize',14);
                elseif strcmp(Resource.DisplayWindow(i).Type,'Verasonics')
                    % Set axes limits base on Resource.DisplayWindow(i).Orientation
                    if Trans.type < 2
                        % DisplayWindow is x,z oriented
                        vmin = Resource.DisplayWindow(i).ReferencePt(3); % vertical minimum is DisplayWindow.ReferencePt(3).
                        vmax = vmin + Resource.DisplayWindow(i).Position(4)*Resource.DisplayWindow(i).pdelta;
                        xmin = Resource.DisplayWindow(i).ReferencePt(1);
                        xmax = xmin + Resource.DisplayWindow(i).Position(3)*Resource.DisplayWindow(i).pdelta;
                   else
                        if ~isfield(Resource.DisplayWindow(i),'Orientation')
                            error(['VSX: Resource.DisplayWindow(%d).Orientation must be defined for a 2D Displaywindow',...
                                   'with Trans.type = 2\n'],i);
                        end
                        switch Resource.DisplayWindow(i).Orientation
                            case 'xz'
                                % DisplayWindow is x,z oriented
                                vmin = Resource.DisplayWindow(i).ReferencePt(3);
                                vmax = vmin + Resource.DisplayWindow(i).Position(4)*Resource.DisplayWindow(i).pdelta;
                                xmin = Resource.DisplayWindow(i).ReferencePt(1);
                                xmax = xmin + Resource.DisplayWindow(i).Position(3)*Resource.DisplayWindow(i).pdelta;
                            case 'yz'  % DisplayWindow is y,z oriented
                                vmin = Resource.DisplayWindow(i).ReferencePt(3);
                                vmax = vmin + Resource.DisplayWindow(i).Position(4)*Resource.DisplayWindow(i).pdelta;
                                xmin = Resource.DisplayWindow(i).ReferencePt(2);
                                xmax = xmin + Resource.DisplayWindow(i).Position(3)*Resource.DisplayWindow(i).pdelta;
                            case 'xy'  % DisplayWindow is x,y oriented
                                vmin = Resource.DisplayWindow(i).ReferencePt(2);
                                vmax = vmin + Resource.DisplayWindow(i).Position(4)*Resource.DisplayWindow(i).pdelta;
                                xmin = Resource.DisplayWindow(i).ReferencePt(1);
                                xmax = xmin + Resource.DisplayWindow(i).Position(3)*Resource.DisplayWindow(i).pdelta;
                        end
                    end

                    % Generate an event that the spatial position has changed.
                    spatialPositionEvent = SpatialPositionEvent(xmin, xmax, vmin, vmax);
                    VantageImageEvent.generateSpatialPositionEvent(Resource.DisplayWindow(i).imageHandle, spatialPositionEvent);
                end
            end
        case 'rcvloop'
            % switch to/from receive data loop mode (Resource.Parameters.simulateMode = 2)
            if rloopButton == 1
                if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                Control(n).Command = 'copyBuffers'; % copyBuffers does a sequenceStop
                runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.
                Resource.Parameters.simulateMode = 2;
                if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                Control(n).Command = 'set&Run';
                Control(n).Parameters = {'Parameters',1,'simulateMode',2,'startEvent',Resource.Parameters.startEvent};
                simButton = 0;
                set(findobj('String','Simulate'),'Value',0);
            else
                if (VDAS==1) % restart sequence
                    Resource.Parameters.simulateMode = 0;
                    if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                    Control(n).Command = 'set&Run';
                    Control(n).Parameters = {'Parameters',1,'simulateMode',0,'startEvent',Resource.Parameters.startEvent};
                    simButton = 0;
                    set(findobj('String','Simulate'),'Value',0);
                else % if no hardware, go back to simulate mode 1
                    Resource.Parameters.simulateMode = 1;
                    if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                    Control(n).Command = 'set&Run';
                    Control(n).Parameters = {'Parameters',1,'simulateMode',1,'startEvent',Resource.Parameters.startEvent};
                    simButton = 1;
                    set(findobj('String','Simulate'),'Value',1);
                end
            end
        case 'simulate'
            % switch mode to/from simulate.
            if (simButton == 1)
                Resource.Parameters.simulateMode = 1;
                if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                Control(n).Command = 'set&Run';
                Control(n).Parameters = {'Parameters',1,'simulateMode',1,'startEvent',Resource.Parameters.startEvent};
                if rloopButton == 1
                    rloopButton = 0;
                    set(findobj('String','Rcv Data Loop'),'Value',0);
                end
            else
                if (VDAS==1)
                    Resource.Parameters.simulateMode = 0;
                    if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
                    Control(n).Command = 'set&Run';
                    Control(n).Parameters = {'Parameters',1,'simulateMode',0,'startEvent',Resource.Parameters.startEvent};
                else
                    simButton = 1; % if no hardware, change state back to simulate.
                    set(findobj('String','Simulate'),'Value',1);
                    if Resource.Parameters.verbose>1
                        fprintf('VSX: Entering simulate mode since hardware is not present.\n');
                    end
                end
            end
    end

    action = 'none';  % this prevents continued execution of same action.
    % if ~isempty(Control(n).Command), disp(Control), drawnow, end % TEST: Uncomment to print out Control changes.

    %% Call VsUpdate and TXEventCheck before calling runAcq
    % Make a copy of Control, and clear it so any new commands will
    % be processed next time.
    VSX_Control = Control;
    Control = struct('Command', [], 'Parameters', []);
    if ~isempty(VSX_Control) && ~isempty(VSX_Control(1).Command)
        % There is a command to be processed.
        % Create flags to track which update calls are needed, so they
        % will be done in the correct order and only once.
        updateTW = false;
        updateTX = false;
        updateReceive = false;
        updateSeqControl= false;
        checkTXLimits = false;
        for ctlnum = 1:length(VSX_Control)
            if strcmp(VSX_Control(ctlnum).Command, 'update&Run')
                for parnum = 1:length(VSX_Control(ctlnum).Parameters)
                    switch VSX_Control(ctlnum).Parameters{parnum}
                        case 'TW'
                            updateTW = true;
                            updateTX = true;
                            checkTXLimits = true;
                        case 'TX'
                            updateTX = true;
                            checkTXLimits = true;
                        case 'Receive'
                            updateReceive = true;
                            checkTXLimits = true;
                        case 'SeqControl'
                            updateSeqControl = true;
                            checkTXLimits = true;
                        case {'Event', 'Parameters', 'Trans'}
                            checkTXLimits = true;
                    end
                end
            elseif strcmp(VSX_Control(ctlnum).Command, 'set&Run')
                for parnum = 3:4:length(VSX_Control(ctlnum).Parameters)
                    if strcmp(VSX_Control(ctlnum).Parameters(parnum), 'startEvent')
                        checkTXLimits = true;
                    end
                end
            end
        end
        % Now process the updates that have been identified.
        if updateTW; updateh('TW'); end
        if updateTX; updateh('TX'); end
        if updateReceive; updateh('Receive'); end
        if updateSeqControl; updateh('SeqControl'); end
        if VDASupdates && checkTXLimits
            % call TXEventCheck to process the updated structures in the
            % matlab workspace
            TXEventCheckh();
            % Add an update&Run command for the TPC structure, so runAcq
            % will read in and process the updated values
            n=length(VSX_Control)+1;
            VSX_Control(n).Command = 'update&Run';
            VSX_Control(n).Parameters = {'TPC'};
        end
        clear updateTW updateTX updateReceive updateSeqControl checkTXLimits ctlnum parnum
    end

    % Perform low pass filter on sequencePeriod times.
    sequencePeriod = sequencePeriod * 0.8 + tElapsed * 0.2;
    % Call runAcq and time its execution.
    tStartRunAcq = tic;
    % NOTE:  If runAcq() has an error, try/catch below will allow closing
    % HW before reporting error and exiting.
    try
        runAcq(VSX_Control);
    catch msg
        disp(' ');
        fprintf(2, 'RUNACQ ERROR:  Closing HW before displaying error and exiting.\n\n');
        if (VDAS == 1)
            % HIFU: before closing HW, check for use of external power supply and disable it.
            if TPC(5).inUse == 2
                % save hv setting so it will be in the workspace after we
                % quit (the following commands will set it back to 1.6)
                lasthv5 = TPC(5).hv;
                % External supply was in use so now we have to disable its output and set voltage and
                % current back to minimum levels of 1.6 V., 2 A.
                [~, ~] = extPwrCtrl('CLOSE',1.6); % disable and close
                TPC(5).hv = lasthv5;
            end
            % Close hardware.
            [~] = hardwareClose();
        end
        rethrow(msg)
    end
    tElapsed = toc(tStartRunAcq);
    if freeze==1    % if freeze got set by runAcq, set state of togglebutton.
        set(findobj('String','Freeze'),'Value',1);
       
    end
    exit=1;
end  % end of the VSX runtime while loop

% Copy buffers back to Matlab space. copyBuffers also stops the hardware sequencer.
Control(1).Command = 'copyBuffers';
runAcq(Control); % NOTE: If runAcq() has an error, it reports it then exits MATLAB.
if (VDAS == 1)
    % HIFU: before closing HW, check for use of external power supply and disable it.
    if TPC(5).inUse == 2
        % save hv setting so it will be in the workspace after we
        % quit (the following commands will set it back to 1.6)
        lasthv5 = TPC(5).hv;
        % External supply was in use so now we have to disable its output and set voltage and
        % current back to minimum levels of 1.6 V., 2 A.
        [~, ~] = extPwrCtrl('CLOSE',1.6); % disable and close
        TPC(5).hv = lasthv5;
    end
    % Close hardware.
    Result = hardwareClose();
end
if Resource.Parameters.verbose>1
    % Print out frame rate estimate.
    if ~exist('frameRateFactor','var')
        fprintf('Sequence rate = %3.2f\n', 1/sequencePeriod)
    else
        fprintf('Frame rate = %3.2f\n', frameRateFactor/sequencePeriod)
    end
end
% If a LogData record size was specified, convert the LogData file.
% *** Note: The called routine is platform specific.
if (ismac && shouldConvertLogData), convertLogData; end

return
