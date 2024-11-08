function F = cbaload(fname, isstd, aconf, bchan, fconf, prof, disp, ver)

% CBALOAD read serum ab cytometric assay FCS file and initialize struct
%
%    A = CBALOAD(fcs, isstd, aconf, bchan, fconf, prof, disp, ver) reads
%    array, flow and peak configuration files and the specified FCS file
%    and returns a data structure used for subsequent peak finding and
%    reactivity determination.  Input arguments are: full path to FCS file,
%    flag indicating regular sample (0) or standard control sample (1), 
%    array configuration file, bead peak detection channel (and split 
%    channel if used) in the FCS file, flow configuration string, peak
%    detection parameter file, flag indicating whether to display (1) or
%    suppress (0) figures, and software version (passed from parent
%    script). This script is typically not run in isolation but called from
%    parent script cbafcs.
%
%    Usage:
%
%        A = CBALOAD('FCS/myflowfile.fcs', 0, 'conf/myarrayconf.txt',...
%                     'FL4-A', 'conf/myflowconf.txt',...
%                     'conf/myprofile.txt', 1, '1.2');
%
%    If using a split channel, specify after the bead channel separated by
%    a colon (in this case, the array config file should be modified
%    accordingly, as well):
%
%        A = CBALOAD('FCS/myflowfile.fcs', 0, 'conf/myarrayconf2.txt',...
%                     'FL4-A:FL5-A', 'conf/myflowconf.txt',...
%                     'conf/myprofile.txt', 1, '1.2');
%
%    Uses the function "fca_readfcs.m" to read FCS files
%    Copyright (c) 2020, Laszlo Balkay, All rights reserved.
%
% -------------------------------------------------------------------------
%
%    Authors: Alexander F. Rosenberg (afr@uab.edu) and Rodney G. King
%        with John T. Killian, Fen Zhou, Davide Botta, Todd J. Green,
%        Jobaida Akther, M. Emon Hossain, Shihong Qiu, Guang Yang,
%        Troy D. Randall and Frances E. Lund
%
%    University of Alabama at Birmingham
%    Department of Microbiology
%    April 11, 2023
%    Copyright (C) 2024 UAB Research Foundation
%    This software is offered with no guarantees of any kind.
%
%    see: "A high-throughput multiplex array for antigen-specific serology
%    with automated analysis", bioRxiv 2023 April.
%    doi: 10.1101/2023.03.29.534777
%
%    This file (part of the "CBA Toolbox") is free software: you can
%    redistribute it and/or modify it under the terms of the GNU General
%    Public License as published by the Free Software Foundation, version 3
%    of the License.  This file is distributed in the hope that it will be
%    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    General Public License for more details.  You should have received a
%    copy of the GNU General Public License along with this program.  If  
%    not, see https://www.gnu.org/licenses/gpl-3.0.en.html.


    % validate input args
    ty = {'string', 'char'};
    ne = {'nonempty'};
    pa = inputParser;
    addRequired(pa, 'fname', @(x) validatefile(x));
    addRequired(pa, 'isstd', @(x) ismember(x, [0 1]));
    addRequired(pa, 'aconf', @(x) validatefile(x));
    addRequired(pa, 'bchan', @(x) validateattributes(x, ty, ne));
    addRequired(pa, 'fconf', @(x) validateattributes(x, ty, ne));
    addRequired(pa, 'prof',  @(x) validatefile(x));
    addRequired(pa, 'disp',  @(x) ismember(x, [0 1]));
    addRequired(pa, 'ver',   @(x) validateattributes(x, ty, ne));
    parse(pa, fname, isstd, aconf, bchan, fconf, prof, disp, ver);
    
    % break down bead channel and split channel if needed: split:bead
    ch = strsplit(char(bchan), ':');
    bc = ch{1};
    if length(ch) == 2, sc = ch{2}; else, sc = NaN; end
       
    % initialize data structure
    F = struct;
    F.version      = char(ver);
    F.analysisdate = char(datetime);
    F.fcsname      = char(fname);
    F.arrayconf    = char(aconf);
    F.profile      = char(prof);
    F.beadchan     = bc;
    F.splitchan    = sc;
    F.flowconf     = char(fconf);    
    F.standard     = isstd;
    F.display      = disp;
        
    % read bead array config (clear spaces around =) - values are strings
    arows = regexprep(readcell(aconf, 'delimiter', '\t'), '\s*=\s*', '=');
    arows = reshape(strsplit(strjoin(arows, '='), '='), 2, size(arows, 1))';
    A = struct;
    for i = 1:length(arows), A.(arows{i, 1}) = arows{i, 2}; end

    % extract ag orders for each bead size (and split bin if specified)
    fa = fields(A);
    sz = table;
    sz.field = string(fa(~cellfun('isempty', regexp(fa, 'bead_size_'))));
    sz.size = str2num(char(regexprep(sz.field, 'bead_size_', '')));
    sz.label = cell(1, length(sz.size))';
    for i = 1:length(sz.size)
        sz.label(i) = cellstr(A.(sz.field{i}));
    end
    sz.label = string(sz.label);
    F.beadsize = sz;
    if isstd == 0, ty = 'peak'; else, ty = 'ctrl'; end
    fb = ['bead_' ty '_order_'];
    beadbin = fa(~cellfun('isempty', regexp(fa, fb)));
    bc = nan(length(beadbin), 4);
    for i = 1:length(beadbin)
        suffix = regexprep(beadbin{i}, fb, '');
        c = str2num(cell2mat(strsplit(suffix, '_')'));
        ag = strsplit(A.(beadbin{i}), ',');
        bc(i, 1) = c(1);
        bc(i, 2) = 1; if length(c) == 2, bc(i, 2) = c(2); end
        bc(i, 3) = length(ag);
        bc(i, 4) = bc(i, 3) - length(find(strcmpi(ag, 'null')));
        F.(['peakorder' regexprep(suffix, '_', 'p')]) = A.(beadbin{i});
    end
    F.beadconf  = array2table(bc, 'VariableNames', ...
        {'size', 'bin', 'npeaks', 'used'});
    F.nbeadsize = length(unique(F.beadconf.size));
    F.nbeadbin  = size(F.beadconf, 1);
    
    % read peak profile (clear spaces around =) - values are numeric
    p = regexprep(readcell(prof, 'delimiter', '\t'), '\s*=\s*', '=');    
    p = reshape(strsplit(strjoin(p, '='), '='), 2, size(p, 1))';
    F.settings = struct;
    for i = 1:length(p), F.settings.(p{i, 1}) = str2num(p{i, 2}); end %#ok<*ST2NM>
    
    % miscellaneous settings (rows, three bead sizes) - columns for splits
    F.pref.colororder = {...
        [1.0 0.0 0.0], [1.0 0.5 0.5];...
        [0.0 0.0 1.0], [0.5 0.5 1.0];...
        [0.0 0.6 0.0], [0.3 0.8 0.3]};
        
    % required parameters: col1 = parameters, col2 = defaults (UNUSED)
    key = {...
        'beadread_asinh_scale_factor',  3000;...        
        'beadsize_xfsc_range',          [3.5 7.5];...
        'beadsize_yssc_range',          [3.5 8.5];...
        'beadsize_oval_gate_scale',     1.2;...
        'beadsize_bin_size',            .01;...
        'beadsize_loess_smooth',        .04;...
        'beadsize_min_peak_height',     .1;...
        'beadsize_min_peak_prominence', .1;...
        'beadsing_dbclus_epsilon',      .02;...
        'beadsing_dbclus_min_points',   50;...
        'beadpeak_range',               [-1 9];...
        'beadpeak_bin_size',            .075;...
        'beadpeak_loess_smooth',        .04;...
        'beadpeak_min_peak_height',     .04;...
        'beadpeak_min_peak_prominence', .04;...
        'beadreac_range',               [-1 7];...
        'beadreac_bin_size',            .05;...
        'beadreac_loess_smooth',        .1};

    % make sure all parameters are in imported profile
    d = setdiff(key(:, 1), fieldnames(F.settings));
    if ~isempty(d)
        error(['parameters not in profile: ' strjoin(d, ', ')]);
    end
        
    % parse isotype channels
    isodef = strsplit(fconf, ',');

    % read fcs file
    [d, h] = fca_readfcs(fname);
    F.header = h;
    F.chan = {h.par.name};
    F.raw = d;    
    
    % asinh-transformed channels of interest
    sf = F.settings.beadread_asinh_scale_factor;
    F.event = table;   
    F.event.SSC_H = asinh(d(:, F.chan == "SSC-H") / sf);
    F.event.FSC_H = asinh(d(:, F.chan == "FSC-H") / sf);
    F.event.SSC_A = asinh(d(:, F.chan == "SSC-A") / sf);
    F.event.FSC_A = asinh(d(:, F.chan == "FSC-A") / sf);

    % error checking for FCS flow channel specifications
    if ~ismember(F.beadchan, F.chan)
        error(['bead channel ' F.beadchan ' not found in FCS header']);    
    elseif ~isempty(setdiff(regexprep(isodef, '.*:', ''), F.chan))
        error('one or more isotype flow channel not in FCS header');    
    elseif ~isnan(F.splitchan) & ~ismember(F.splitchan, F.chan)
        error(['split channel ' F.splitchan ' not in FCS header']);
    end
    
    % assign & transform flow config channels: bead channel
    F.event.bead = asinh(d(:, F.chan == string(F.beadchan)) / sf);

    % assign & transform flow config channels: split channel
    if ~isnan(F.splitchan)
        F.event.split = asinh(d(:, F.chan == string(F.splitchan)) / sf);
    end
    
    % assign & transform flow config channels: secondaries channels
    for i = 1:length(isodef)
        a = strsplit(isodef{i}, ':');
        F.event.(a{1}) = asinh(d(:, F.chan == string(a{2})) / sf);
    end
    
return
    
    
