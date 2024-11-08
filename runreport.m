function D = runreport(A, expt, varargin)

% RUNREPORT generate batch report with optional concentration calculation
%
%    RUNREPORT(X, ex) takes as input a structure array created by RUNEXPT
%    as the first argument and an experiment name as the second 
%    argument and returns a structure array of tables (one per flow 
%    config) of raw MFI data, arcsinh-transformed MFI data and optionally,
%    back-calculated concentrations (arbitrary concentration units).
%    Additionally, output tables are written to TXT file(s) - one per flow
%    configuration.
%    
%    RUNREPORT(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%    parameter/value pair values:
%        'stdcurve'     specify a standard curve object (as created by the
%                       function RUNSTD) - this will contain the four
%                       curve fit parameters for each isotype, bead size
%                       and low configuration
%        'outdir'       specify a directory for the TXT file to be written
%                       to - default is the directory from which this
%                       function is run
%        'txt'          flag to determine whether to write tables to output
%                       text files.  values can be "yes" or "no" - default
%                       if unspecified is "yes"
%
%    Usage:
%
%        Y = RUNREPORT(A, 'myexpt');
%
%        Y = RUNREPORT(A, 'myexpt', 'outdir', '/Users/afr/output/');
%
%        Y = RUNREPORT(A, 'myexpt', 'txt', 'no');
%
%        Y = RUNREPORT(A, 'myexpt', 'stdcurve', Q, 'outdir', './results/');
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
%    November 7, 2024
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


    % input parameter-value error checking
    pa = inputParser;
    validateattributes(expt, {'string', 'char'}, {'nonempty'});
    addParameter(pa, 'outdir',  './',...
        @(x) validateattributes(x, {'string', 'char'}, {'nonempty'}));
    addParameter(pa, 'txt', 'yes',...
        @(x) validateattributes(x, {'string', 'char'}, {'nonempty'}));
    addParameter(pa, 'stdcurve', 0);
    parse(pa, varargin{:});    
    outdir = char(pa.Results.outdir);
    S = pa.Results.stdcurve;
    expt = string(expt);
    if ~isequal(outdir(end), '/'), outdir = [outdir '/']; end
    outdir = string(outdir);
    if lower(string(pa.Results.txt)) == "yes"
        outfile = 1;
    else
        outfile = 0;
    end 

    % get non-control samples
    M = A(~[A.iscontrol]);
    n = length(M);

    % build output table from non-control data
    Z = table;
    Z.experiment = repmat(expt, n, 1);
    Z.donor = string({M.donor}');
    Z.day = [M.days]';
    Z.fcs = string({M.fcs}');
    Z.plate = string({M.plate}');
    Z.well = string({M.well}');
    Z.dilution = [M.dilution]';
    % Z.iscontrol = [M.iscontrol]';
    % Z.concentration = [M.concentration]';
    Z.flowconfig = repmat("", n, 1);

    % add flowconfigs (per sample)
    for i = 1:n, Z.flowconfig(i) = string(M(i).results.flowconf); end
    
    % feature base names (use based on first sample, assume all the same)
    aconf = num2str(M(1).results.out.size) + "_" + ...
        num2str(M(1).results.out.peak) + "_" + ...
        regexprep(M(1).results.out.name, ' ', '_');
    
    % add # events per feature, per sample
    events = nan(n, length(aconf));
    for i = 1:n, events(i, :) = M(i).results.out.events'; end
    E = array2table(events, 'VariableNames', "events_" + aconf);
    Z = [Z E];

    % are dilutions and flowconfigs always linked?
    ufc = unique(Z.flowconfig);

    % separate table by flowconfig
    D = struct;
    for i = 1:length(ufc)
        D(i).flowconf = ufc(i);
        D(i).idx      = find(Z.flowconfig == ufc(i));
        D(i).M        = M(D(i).idx);
        D(i).table    = Z(D(i).idx, :);
    end

    % are FCS files guaranteed to be unique ???
    
    % assume asinh scale factor same for all samps (if created by runexpt)
    sf = M(1).results.settings.beadread_asinh_scale_factor;
    
    % step through flowconfigs
    nfeat = length(aconf);
    for i = 1:length(ufc)

        isotypes = regexprep(split(ufc(i), ','), ':.*$', '')';
        niso = length(isotypes);
        nsamp = length(D(i).idx);
        h = reshape(repmat(isotypes, length(aconf), 1) + "_" + repmat(aconf, 1, niso), length(aconf) * niso, 1)';
        F = nan(nsamp, nfeat * niso);
        for j = 1:nsamp
            temp = nan(nfeat, niso);
            for k = 1:niso, temp(:, k) = D(i).M(j).results.out.(isotypes(k)); end
            F(j, :) = temp(:)';
        end
        F2 = array2table(F, 'VariableNames', "asinhMFI_" + h);
        F3 = array2table(sinh(F) * sf, 'VariableNames', "rawMFI_" + h);
        D(i).table = [D(i).table F3 F2];

        % if a standard curve structure was specified
        if isstruct(S)
    
            % make sure flowconfigs match that of standard curve object
            ufcq = [S.flow]';
            if ~isequal(sort(ufc), sort(ufcq))
                error('flowconfig in data and stdcurve object not matched');
            end

            % dilution
            df = D(i).table.dilution;

            % standard curve fits
            K = S(ufc(i) == ufcq).fitparam;            

            % get isotypes and bead sizes for each column of mfi data
            sz = M(1).results.out.size;
            H = table;
            H.iso = reshape(repmat(isotypes, nfeat, 1), niso * nfeat, 1);
            H.sz = repmat(sz, niso, 1);

            % make blank concentration table
            C = nan(size(F));
            
            % iterate through isotypes and bead sizes, compute concentrations
            for iiso = 1:length(isotypes)
                for isz = 1:length(unique(K.beadsize))
                    iparam = find(K.isotype == isotypes(iiso) & K.beadsize == isz);
                    pa = K.A(iparam);
                    pb = K.B(iparam);
                    pc = K.C(iparam);
                    pd = K.D(iparam);
                    icol = find(H.iso == isotypes(iiso) & H.sz == isz);
                    m = F(:, icol);
                    temp0 = (10.^(pc * ((((pa - pd) ./ (m - pd)) - 1) .^ (1 / pb))));
                    correction = repmat((df / 1000), 1, length(icol));
                    temp = temp0 .* correction;  % ug/ml
                    temp(m < pa) = -1e99;
                    temp(m > pd) = 1e99;
                    C(:, icol) = temp;            
                end
            end
            C2 = array2table(C, 'VariableNames', "arbConc_" + h);
            D(i).table = [D(i).table C2];

        end

        if outfile == 1

            % file name
            dname = string(datetime('today', 'format', 'yyyyMMdd'));
            aname = strsplit(A(1).results.arrayconf, '/');
            aname = regexprep(aname{end}, '\.txt', '');
            fname = regexprep(char(ufc(i)), {':', ','}, {'_', '_'});
            oname = outdir + expt + "_" + dname + "_" + aname + "_" + fname + ".txt";
    
            % write to output file
            writetable(D(i).table, oname, 'FileType', 'text', 'Delimiter', '\t');

        end
    end

    % clean up return value
    D = rmfield(D, {'flowconf', 'idx', 'M'});

return




    
