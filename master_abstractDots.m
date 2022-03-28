function master_abstractDots(abstractData)

% This functions generates plots from Shushruth et al. 2022 (Curr. Biology
% The relevant data file is present in the \Data folder
% It is also available on Mendeley Data at DOI:
% Please add the folder \Modules to the path for supporting functions

% Data in the field "annie" is from Monkey-AN and 
% data in the field smith is from Monkey-SM

% The data from each monkey is contained in the following fields
% .annie.eMat: behavioral data from monkey-AN
% .annie.respMat: neural data from monkey-AN
% .annie.ME: motion energy data from monkey-AN
% .smith.eMat: behavioral data from monkey-SM
% .smith.eMatRev: behavioral data from monkey-SM with target colors 
%           reversed for cells [10 13 15 19 23 24] (see STAR Methods)
% .smith.respMat: neural data from monkey-SM
% .smith.ME: motion energy data from monkey-SM

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries and set up
dataOut = struct();
E = eventCodes;
annie = struct(); % pre allot
smith = struct(); % pre allot

anMat = abstractData.annie.eMat;
smMat = abstractData.smith.eMat;
smMatRev = abstractData.smith.eMatRev;

anRespMat = abstractData.annie.respMat;
smRespMat = abstractData.smith.respMat;

anME = abstractData.annie.ME;
smME = abstractData.smith.ME;

%% Figure 2 A & B
annie.behavior = abstractDots_psychometrics(anMat,1,'Fig. 2A (AN)');
smith.behavior = abstractDots_psychometrics(smMat,2,'Fig. 2B (SM)');
