%% Load data

load("IOdownload.mat");

%% 1a) The in-degree and out-degree centrality

% Input matrices
WSwe = io.swe2000;
WIdn = io.idn2000;

% The three most central sectors of Sweden
[in_central_name_swe, out_central_name_swe] = In_out_degree_centrality(WSwe, name);

disp('In-degree centrality of Sweden during year 2000:');
disp(in_central_name_swe);

disp('Out-degree centrality of Sweden during year 2000:');
disp(out_central_name_swe);

% The three most central sectors of Indonesia
[in_central_name_idn, out_central_name_idn] = In_out_degree_centrality(WIdn, name);

disp('In-degree centrality of Indonesia during year 2000:');
disp(in_central_name_idn);

disp('Out-degree centrality of Indonesia during year 2000:');
disp(out_central_name_idn);

%% 1b) The eigenvector centrality on the largest connected component

% The three most central sectors of Sweden
sectors_swe = Eigenvector_centrality(WSwe,name);

disp('Eigenvector centrality of Sweden during year 2000: ')
disp(sectors_swe)

% The three most central sectors of Indonesia
sectors_idn = Eigenvector_centrality(WIdn,name);

disp('Eigenvector centrality of Indonesia during year 2000: ')
disp(sectors_idn)

%% 1c) The Katz centrality

% The three most central sectors of Sweden
[sectors_swe, sectors_swe_wrtr] = Katz_centrality(WSwe,name);

disp('Katz centrality of Sweden during year 2000: ')
disp(sectors_swe)

disp('Katz centrality of Sweden during year 2000, Wholesale & retail trade; repairs: ')
disp(sectors_swe_wrtr)

% The three most central sectors of Indonesia
[sectors_idn, sectors_idn_wrtr] = Katz_centrality(WIdn,name);

disp('Katz centrality of Indonesia during year 2000: ')
disp(sectors_idn)

disp('Katz centrality of Indonesia during year 2000, Wholesale & retail trade; repairs: ')
disp(sectors_idn_wrtr)