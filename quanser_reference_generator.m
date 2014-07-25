%QUANSER_REFERENCE_GENERATOR Generate a reference path using the Simulink
%model 'reference_generator.mdl'
%% Define model specfic data
N = 600; % samples
nx = 6; % states
nu = 2; % inputs
uech = [1.8; 1.8]; % The stabilizing command
matfilename = 'references/ref1.mat';
%% Run the reference generation Simulink model
%This creates a timeseries variable, simout, with the reference paths.
%The first colums is the elevation reference and second is the pitch
%reference
sim('reference_generator');
data = simout.Data(1:N, :)';
XREF = zeros(nx, N);
UREF = repmat(uech, [1 N]);
XREF(1,:) = data(1,:); % import the elevation reference
XREF(3,:) = data(2,:); % import the pitch reference
%% Write file and clean-up
save(matfilename,'XREF','UREF','-v7');
clear simout data