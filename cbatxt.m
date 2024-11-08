function cbatxt(F, outdir)

% CBATXT prepare output info for TXT file
%    CBATXT(A, out) generates a text file report based on the output of
%    cbafcs. The first argument is a complete data structure, and the 
%    second argument is the output directory.
%
%    Usage:
%
%        CBATXT(A, '~/user/');
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

    
    % get header
    H = cbaheader(F);

    % complete header
    htxt = [H.part1; H.part2; H.part3; H.part4; H.part5];
    htxt2 = cell(size(htxt, 1), 1);
    for i = 1:size(htxt, 1)
        htxt2(i) = cellstr([htxt{i, 1} '=' htxt{i, 2}]);
    end

    % reactivity table (add a bin column if a split channel was used)
    [nrow, ncol] = size(F.out);   
    fn = F.out.Properties.VariableNames;
    header = cellstr(strjoin(fn, '\t'));
    str = '%d\t%d\t%f\t%f\t%s\t%d\t';
    if ~isnan(F.splitchan), str = ['%d\t' str]; end
    str = [str repmat('%f\t', 1, ncol - find(strcmp(fn, 'events')))];
    str = str(1:(end - 2));
    tab = cell(nrow, 1);
    for i = 1:nrow
        x = table2cell(F.out(i,  :));
        tab(i, 1) = cellstr(sprintf(str, x{:}));
    end    
    lines = [htxt2; header; tab];        

    % base file name
    fcs_name = F.fcsname(max(strfind(['/' F.fcsname], '/')):end);
    fcs_base = fcs_name(1:(max(strfind(fcs_name, '.')) - 1));
    
    % write to file
    fid = fopen([outdir fcs_base '_OUTPUT.txt'], 'w');
    for i = 1:length(lines), fprintf(fid, '%s\n', lines{i}); end
    fclose(fid);

