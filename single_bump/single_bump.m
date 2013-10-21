% [r, u] = single_bump()
%
% Runs code for a bump attractor on a linear space, following homework 2 of Neu
% 437, Spring 2011. All parameters for the simulation have built-in default
% values, so the code runs even when called with no arguments, e.g.:
%
%  >> single_bump
%
% PASSING IN OPTIONAL PARAMETERS: You can change parameters from their
% default values either through the graphical user interface (GUI)
% controls provided in the simulation's window, or at the time of calling
% single_bump.m by using a (parameter_name, parameter_value) pair.  For
% example, to specify that Ncells should be 150 and the simulation timestep
% should be 0.05, you pass in two pairs and call:
%
%  >> single_bump('Ncells', 150, 'dt' 0.05)
%
% Parameters not specified in your list of (name, value) pairs will be kept
% at their default values. The order of (parameter_name, parameter_value)
% pairs doesn't matter, you can put them in any order when calling
% single_bump.m
%
% GRAPHICAL MODE VS. NON-GRAPHICAL MODE: If the optional parameter
% 'do_plot' is passed as 0, then no graphics or GUI are generated. In that
% case, the code runs much faster. The return values (see below) will be
% the r vectors u at the end of the simulation.
%
% KEYSTROKES: In graphical mode, click on the main axis and then press
% 'p' to produce a positive input jolt to the network at the mouse's
% location. Press 'd' to produce a negative jolt. (At least on a Mac, you
% can keep the key pressed and move the mouse to produce a rapid sequence
% of jolts.) Press 'q' to end the graphical simulation.
%
% RETURNS:
% --------
%
%   r      A vector, Ncells long, that represents the firing rate of each
%          neuron at the end of the simulation.
%
%   u      A vector, Ncells long, that represents the membrane potential of
%          each neuron at the end of the simulation.
%

% written by Carlos Brody for Neu 501, Fall 2009. Modified by


function [rr, u] = single_bump(varargin)


t = 0;        % current time


%% Some initial values

% Parameters of the sigmoid activation function of each cell:
cell_thresh = 0.5;  % r = 0.5 + 0.5*tanh((u - cell_thresh)/cell_sigma)
cell_sigma  = 0.3;

pairs = { ...
   'dt'           0.1   ; ...  % timestep size
   'T'          10000   ; ...  % time after which simulation stops
   'Ncells'       250   ; ...  % number of cells in the simulation
   'wE'           1.4   ; ...  % local excitatory weight strength
   'sigmaE'        4    ; ...  % width of excitatory neighborhood
   'gI'           0.4   ; ...  % global inhibition: from every cell to every cell
   'global_drive' 2.1   ; ...  % an external constant drive to every cell
   'sigma_noise'  0.4   ; ...  % magnitude of noise in each cell
   'leak_noise'    0    ; ...  % Static randomness, across cells, in their "leak membrane conductance"
   'wDA'           0    ; ...  % Strength of external gaussian drive
   'wDB'           0    ; ...  % Strength of external gaussian drive; B track
   'sigmaD'        4    ; ...  % External gaussian drive has same width as default exc conn width
   'muD'          0.5   ; ...  % Drive starts at midpoint of screen.
   'muDB'          0.5  ; ...  % Drive B starts at midpoint of screen. ... permuted?
   'dVA'           0    ; ...  % Drive speed, in neuron positions per unit time.
   'dVB'           0    ; ...  % Drive speed, in neuron positions per unit time; track B
   'initialU'     []    ; ...  % Optional initial conditions for u
   'do_plot'       1    ; ...  % If this is 0, no figure or plot is generated (runs faster).
   'B_order'       0    ; ...  % permuation for field B
   'use_field_B'   1   ; ...  % set to 0 to use only field A
}; parseargs(varargin, pairs);

if use_field_B
    use_field_B = 1;
end
muD = muD*Ncells;
muDB = muDB*Ncells;

leak = ones(Ncells,1) + leak_noise*randn(Ncells,1); leak(leak<0) = 0;  %#ok<NODEF>

if ~isempty(initialU)  %#ok<NODEF>
   if numel(initialU) ~= Ncells,
      error('If you use initialU it must be a vector Ncells long');
   end;
   if size(initialU,1)==1, initialU = initialU'; end;
end;
   


%% Set up the figure
if do_plot,
   figure(1400); clf;
   pos = get(gcf, 'Position');
   % Make sure figure is of enough width, and tell it what fn to call if user
   % presses key:
   set(gcf, 'Units', 'pixels', 'Position', [pos(1) pos(2) 750 800], ...
      'KeyPressFcn', 'bump_keypress_callback');
   % Use the following line instead of prev one if you want two axes stacked
   % on top of each other:400
   % set(gcf, 'Units', 'pixels', 'Position', [pos(1) pos(2) 750 600], ...
   %  'KeyPressFcn', 'bump_keypress_callback');
   
   f2 = figure(2800); clf;
   figure(1400);
   drawnow;
   pos = get(gcf, 'Position');
   
   %% Set up the GUI controllable elements
   sigma_noise_control = add_control('sigma_noise', sigma_noise, 0.05, 0.05, pos(3)-200, 20);
   gI_control          = add_control('gI',          gI,          0.05, 0.05, pos(3)-200, 40);
   wE_control          = add_control('wE',          wE,          0.05, 0.05, pos(3)-200, 60);
   global_drive_control= add_control('global_drive',global_drive,0.25, 0.25, pos(3)-200, 100);
   leak_noise_control  = add_control('leak_noise',  leak_noise,  0.05, 0.05, pos(3)-200, 140);
   
   muD_control         = add_control('drive_position', muD,      0.05, 0.05, pos(3)-200, 260);
   dV_control          = add_control('drive_speed',    dVA,       0.05, 0.05, pos(3)-200, 280);
   wD_control          = add_control('drive_strength', wDA,       0.05, 0.05, pos(3)-200, 300);

   muDB_control         = add_control('drive_pos B', muDB,      0.05, 0.05, pos(3)-200, 180);
   dVB_control          = add_control('drive_spd B',    dVB,     0.05, 0.05, pos(3)-200, 200);
   wDB_control          = add_control('drive_str B', wDB,       0.05, 0.05, pos(3)-200, 220);
   
end;

%%  Set up cells and connections between them
% u is the "membrane potential"; r is the "firing rate"
if isempty(initialU), u = zeros(Ncells,1); 
else                  u = initialU;
end;
r = zeros(Ncells,1); 
eWeights = zeros(Ncells, Ncells);  % excitatory weights

% Set up weights:
if B_order==0
    B_order = randperm(Ncells);
end
P = permMat(B_order);

for i=1:Ncells,
  for j=1:Ncells,
    if i~=j,  % no self-connections
      deltaA = abs(i-j);
      deltaB = abs(B_order(i) - B_order(j));
      eWeights(i,j) = eWeights(i,j)+(...
          exp(-deltaA.^2/(2*sigmaE^2)) + ...
          exp(-deltaB.^2/(2*sigmaE^2)) * use_field_B...
          )/sigmaE;
    end;
  end;
end;


%% Now to plotting and looping
% set up plotting axes:
if do_plot
    figs = [figure(1400), figure(2800)];
end
orders = {1:Ncells, B_order};
wDs = [wDA, wDB];
mus = [muD, muDB];
handles1_ = [];
handles2_ = [];
axes1_ = [];
axes2_ = [];
titles = {'field A: ', 'field B: '};
if do_plot,
    for i=1:2
    order = orders{i};
    w = wDs(i);
    mu = mus(i);
    
       fig = figs(i);
       figure(fig)
       axes2_(i) = axes('Position', [0.1 0.15 0.6 0.1]); set(axes2_(i), 'Units', 'pixels');
       external_drive = w*exp(-(order'-mu).^2/(2*sigmaD.^2));
       handles2_(i) = plot(external_drive, 'k.-'); ylim([-0.5 0.5]);
       xlabel('cell #');

       axes1_(i)  = axes('Position', [0.1 0.3 0.6 0.6]); set(axes1_(i), 'Units', 'pixels');
       handles1_(i) = plot(r, '.-'); set(handles1_(i), 'Color', [0 0.5 0]); ylim([-0.05 1.05]);
       title(titles{i});
       ylabel('activity');
    end;
end

if do_plot
    figure(1400);
end


% Loop forever -- or until user presses 'q'
while t<T,
  t = t + dt;
  
  if do_plot,
     % The following line is to deal with reading keypresses-- no need to modify
     
     % I have no idea what this 1000 handle is, and it's now broken so I'll
     % comment it out. All these global variables and stateful plot handles
     % are really bad form, and Mathworks should feel bad.
     
%      c = get(1000, 'UserData'); if isempty(c), c = ''; end; set(1000, 'UserData', '');
%      % End reading keypress code
%   
%      % c will be either empty if no keypress, or will be the character pressed
%      switch c,
%         % 'p' = Positive Pulse;  'd' = Negative Pulse
%         case {'p' 'd'},
%            [xvalue0, yvalue] = get_pointer_location(ax1);  %#ok<NASGU> % this gets the x position of pointer in axes units
%            for xvalue = xvalue0-1:xvalue0+1,
%               if xvalue >=1 && xvalue <=Ncells, % There are no cells outside the axes
%                  switch c,
%                     case 'p', u(xvalue)=1;  % user asked for positive jolt
%                     case 'd', u(xvalue)=-1; % user asked for negative jolt
%                  end;
%               end;
%            end;
%            % 'q' = quit
%         case 'q',
%            break;  % if q, get out of the WHILE loop
%         otherwise,
%      end;

     % Get the current GUI values of the following variables
     sigma_noise = get_control(sigma_noise_control);
     gI          = get_control(gI_control);
     wE          = get_control(wE_control);
     wDA          = get_control(wD_control);
     dVA          = get_control(dV_control);
     muD         = get_control(muD_control);

     wDB         = get_control(wDB_control);
     dVB         = get_control(dVB_control);
     muDB        = get_control(muDB_control);
     
     global_drive= get_control(global_drive_control);
     
     new_leak_noise = get_control(leak_noise_control);
     if new_leak_noise ~= leak_noise,  % leak_noise doesn't change from timestep to timestep
        leak_noise = new_leak_noise;    % Only if the user changed it do we calculate it anew.
        leak = ones(Ncells,1) + leak_noise*randn(Ncells,1); leak(leak<0) = 0;
     end;
  end;
  
  
  % ---- Dynamics of driving ---
  muD = muD + dVA*dt;
  muDB = muDB + dVB*dt;
  if muD > Ncells-2*sigmaE, dVA=-dVA; if do_plot, set_control(dV_control, dVA); end; end;
  if muD < 2*sigmaE,        dVA=-dVA; if do_plot, set_control(dV_control, dVA); end; end;

  if muDB > Ncells-2*sigmaE, dVB=-dVB; if do_plot, set_control(dVB_control, dVB); end; end;
  if muDB < 2*sigmaE,        dVB=-dVB; if do_plot, set_control(dVB_control, dVB); end; end;
  
  
  if do_plot, set_control(muD_control, muD); end;
  if do_plot, set_control(muDB_control, muDB); end;
  external_drive = wDA*exp(-((1:Ncells)'-muD).^2/(2*sigmaD.^2)) + ...
                   wDB*exp(-((B_order)'-muDB).^2/(2*sigmaD.^2)) * use_field_B;
               
  external_driveB = wDA*exp(-((B_order)'-muD).^2/(2*sigmaD.^2)) + ...
                    wDB*exp(-((1:Ncells)'-muDB).^2/(2*sigmaD.^2)) * use_field_B;
  % ---- End dynamics of driving ---

  
  % ---- Main dynamical equations ----
  dudt = -leak.*u - gI*sum(r) + wE*eWeights*r - sigma_noise*randn(size(u))*sqrt(dt) + ...
         global_drive + external_drive + use_field_B * external_driveB;
  u    = u + dt*dudt;
  r    = 0.5+0.5*tanh((u-cell_thresh)/cell_sigma);
  rB   = 0.5+0.5*tanh(((P*u)'-cell_thresh)/cell_sigma) * use_field_B;
  t = t+dt;
  % ---- End main dynamical equations ----
  
  
  % replot
  external_drives = {external_drive, external_driveB};
  rs = {r, rB};
  for i=1:2
    if do_plot,
        t1 = title(axes1_(i), strcat(...
                   titles{i}, ...
                   sprintf('t=%.2f', t) ...
               ));
        set(handles1_(i), 'YData', rs{i});
        set(handles2_(i), 'YData', external_drives{i});
        drawnow; 
      end;
  end
end;  % and loop

if nargout > 0, rr = r; end;