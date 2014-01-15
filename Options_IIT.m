%% Old options:

op_single = network.options(4);     % just needed for console output
op_context = network.options(6);    % just needed in make_title
op_min = network.options(9);        % just needed in make_title
op_console = options(10);
op_big_phi = network.options(11);
op_normalize = options(13);         % big phi
op_normalize = network.options(14); % small phi

op_complex = network.options(15);

op_small_phi = network.options(16);
op_big_phi_dist = options(17);      
op_average = network.options(18);
op_parallel = in_options(19);


%% new options
op_parallel = in_options(1);
op_average = network.options(2);
op_complex = network.options(3);
op_small_phi = network.options(4);
op_big_phi = network.options(5);

% Not currently normalizing, but option's there for when we might want to
op_normalize = network.options(6); % small phi

op_normalize = options(7);         % big phi

% not necessary:
op_console = options(8);

% Decide whether to use parfor or for loop (which uses precomputed
% matrices) since parfor loop cannot get the network.mat (?)
op_parfor = network.options(9); % used by Animat program

% Whether everything is computed, or we check the connectivity of the
% network beforehand (sometimes people don't have the required
% package/toolbox)
% Also sometimes one might want to compute everything despite the fact that
% it won't form a complex
op_strongconn = network.options(10);

op_extNodes = network.options(11);
op_single = network.options(12);    % just needed for console output == 1