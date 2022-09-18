
% simulate long vs short loops

saveRoot = "T:\2022-02-24_LiangFu\StripeSim_ReelIn3\"; % Update this!
nTrials = 200;
offset = 0;
batchSize = 2;
hidden = true;
GPU_ID = "0"; % "0" or "1"

cmdRuns = cell(0,2);

% general notes: increaseing LEF density strengthens TADs
%    increasing CTCF capture probability strengthens Loops and stripes

% ====  Common Polymer properties ===============
monomers = 1e3; % number of monomers in simulation % monmoer = 1 kb
replicates = 1;  % number of replicated detached chromosomes
trajectoryLength = 5000; % simulation length
density = 0.001;  % monomer density of confining box
% ===== Loop Extrusion parameters =================
% -------- Cohesin --------
% long loops = 330,330.  short loops = 66,66;
LEF_lifetime = 800;%  200;  %   average lifetime before spontanteous fall off. Mirny lab default is 200. could explore this.   
LEF_separation = 500;%  100; % 1e9  ave distance between LEFs (loop extruding factors), measured in monomers. 
loadSites = []; % cohesin load sites, 0 indexed for python
% --------- CTCF ----------


% CTCF 
%  2 creates a stripe pointing right (stops leftward travel)
%  1 creates a stripe pointing left (stops rightward travel);
% 0 points both ways



loadSites =  [205,795]; % for reel model
ctcfSites   =  [200, 800]; % (0 indexed) Position of CTCF sites, recorded in units of monomers
ctcfCapture =   [.999, .999]; % .9 was 'default' probability cohesin is captured rather than walking right past a CTCF site
ctcfRelease =  [.003 .003]; % .003 was 'default', rate constant for cohesin release (smaller numbers hold cohesin for longer)
ctcfDir    =    [2  1]; % 0 is bidirectional, 1 stops right, 2 stops left

% ===== Sticky 'epigenetic' properties ============
monomerTypes = zeros(1,monomers); % initialize
a = []; % % one could load and threshold a ChIP-seq track here to populate this data
monomerTypes(a) = 1; % edit the state assigment of select monomers
interactionMatrix = "[[0,0,0],[0,0,0],[0,0,0]]"; % interaction matrix in python numpy notation (square braces group rows)
% interaction matrix specifies stickiness to self and other types
%    Numbers between 0.1 and 0.4 are typical in Mirny lab simulations
%    Currently configured for up to 3 distinct sticky states.   Diagonal
%    values indicate self-stickiness, off-diagonal indicate other
%    stickiness. Matrix must be symmetric and keep the string notation
%    below.  The monomer states are specified in the "monomerTypes"
% Note, state assignment won't have any effect if the interaction matrix is
% the same for all states. 



% ---- A little house keeping on loop extrusion
% convert load sites to a load probability distribution
%  NOTE: Alterantively you may directly specify a loading probability dist.
%  For example, a ChIP seq track 
loadProb = zeros(1,monomers);
loadProb(loadSites(1)) = .1; % one could instead chose a ChIP-seq track
loadProb(loadSites(2)) = .9; % one could instead chose a ChIP-seq track
loadProb = loadProb + .0001; 
loadProb = loadProb./sum(loadProb);

%======= End of parameter choices ====% 

% location of the python code
polychrom = 'C:\Shared\polychrom-shared\'; 
LEsim = [polychrom,'Simulations\StickyLoopExtSimRpt.py'];


% ====== Run the simulation in python ============
% first we do a little format conversion and the write a .txt file of all
% the parameter choices that will be passed to the python script.
% The code will execute in an external cmd window, freeing up the matlab
% command prompt.  The structure "proc" is returned to matlab and tracks
% the progress of the simulation. For example, proc.HasExited will be false
% (0) while the simulation is running and switch to true (1) when the
% simulation is complete. 
% 
clc
% convert formats. Must be a comma separated string
ctcfSites = num2str(ctcfSites,'% d'); % convert to string (not char array)
ctcfSites = string(['[',regexprep(ctcfSites,' ',','),']']);
for i=1:6
    ctcfSites = regexprep(ctcfSites,',,',',');
end
% ctcfCapture
ctcfCapture = num2str(ctcfCapture,'% .5f'); % convert to string (not char array)
ctcfCapture = string(['[',regexprep(ctcfCapture,' ',','),']'])
% ctcfRelease
ctcfRelease = num2str(ctcfRelease,'% .5f'); % convert to string (not char array)
ctcfRelease = string(['[',regexprep(ctcfRelease,' ',','),']'])
% ctcfDir
ctcfDir = num2str(ctcfDir,'% d'); % convert to string (not char array)
ctcfDir = string(['[',regexprep(ctcfDir,' ',','),']'])
% loadProb
loadProb = num2str(loadProb,'% .8f'); % convert to string (not char array)
loadProb = string(['[',regexprep(loadProb,' ',','),']'])

lp = eval(loadProb);
figure(3); clf; plot(cumsum(lp))
%%
% monomerTypes also needs to be converted to a comma separated string
monomerTypes = string(num2str(monomerTypes)); % convert to string (not char array)
monomerTypes = regexprep(monomerTypes,'  ',',');


simsLaunched = 0;
for t=1:nTrials
    disp(['Running trial ',num2str(t), ' of ',num2str(nTrials)]);
    saveFolder = string([char(saveRoot),'rep',num2str(t+offset,'%03d'),filesep]);

    % write parameter table to pass to the python code
    parsTable = table(monomers,replicates,LEF_lifetime,LEF_separation,ctcfSites,ctcfCapture,ctcfRelease,interactionMatrix,saveFolder,monomerTypes,trajectoryLength,density,loadProb,ctcfDir,GPU_ID);
    parsFile = strcat(saveFolder,'simPars.txt'); % 
    if ~exist(saveFolder)
        mkdir(saveFolder)
    else
        disp('exists, skipping too next');
        continue
    end

    writetable(parsTable,parsFile,'delimiter',';');
    parsPython = char( regexprep(parsFile,'\','/'));
    cmdtxt = ['python "',LEsim,'" "',parsPython,'"'];
    disp(cmdtxt);
    proc = SystemRun(cmdtxt,'Hidden',hidden); % this command actually executes the code    
    simsLaunched = simsLaunched + 1;
    cmdRuns{simsLaunched,1} = proc;
    cmdRuns{simsLaunched,2} = saveFolder;
    
    simsCompleted = sum(cellfun(@(x) x.HasExited, cmdRuns(:,1)));
    simsRunning = simsLaunched - simsCompleted;
    while simsRunning >= batchSize
        pause(10);
        simsCompleted = sum(cellfun(@(x) x.HasExited, cmdRuns(:,1)));
        simsRunning = simsLaunched - simsCompleted;
    end
end
   