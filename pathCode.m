%%
if ispc                                                                     % windows
    addpath(genpath('.\solvers'))
    addpath(genpath('.\matrices')) 
    addpath(genpath('.\referenceSolutions')) 
    addpath(genpath('.\MultiParEig')) 
else                                                                        % not windows
    addpath(genpath('./solvers'))                
    addpath(genpath('./matrices')) 
    addpath(genpath('./referenceSolutions')) 
    addpath(genpath('./MultiParEig')) 
end