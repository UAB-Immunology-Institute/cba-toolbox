function cbapdf(F, outdir)

% CBAPDF prepare output info for PDF file
%    CBAPDF(A, out) generates a text file report based on the output of
%    cbafcs. The first argument is a complete data structure, the 
%    second argument is the output directory.
%
%    Usage:
%
%        CBAPDF(A, '~/user/');
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

    % are 4 or 5 figures made? (i.e. is there a split channel figure)
    usesplit = 0; if ~isnan(F.splitchan), usesplit = 1; end
    
    % base file name
    fcs_name = F.fcsname(max(strfind(['/' F.fcsname], '/')):end);
    fcs_base = fcs_name(1:(max(strfind(fcs_name, '.')) - 1));    

    % import classes for report generator
    import mlreportgen.report.*
    import mlreportgen.dom.*
    
    % remove full paths from fcs file, array config, profile
    k(1) = find(strcmp(H.part1(:, 1), 'FCS_file'));
    k(2) = find(strcmp(H.part1(:, 1), 'array_config'));
    k(3) = find(strcmp(H.part1(:, 1), 'profile'));
    for i = 1:3
        H.part1(k(i), 2) = regexprep(H.part1(k(i), 2), '.*\/', '');
    end   
                    
    % make Table object for main info
    M = Table(H.part1);
    M.Style = [M.Style, {...
        Width('100%'),...
        BackgroundColor('antiquewhite')}];
    M.TableEntriesStyle = {...
        VAlign('middle'),...
        HAlign('left'),...
        InnerMargin('4pt', '4pt', '4pt', '4pt'),...
        FontSize('10pt'),...
        FontFamily('Arial')};    
                        
    % make Table object for parameters
    PR = Table([H.part2; H.part3; H.part4; {'', ''; '', ''}; H.part5]);
    PR.Style = [PR.Style, {...
        Width('100%')}];
    PR.TableEntriesStyle = {...
        VAlign('middle'),...
        HAlign('left'),...
        InnerMargin('0pt', '1pt', '1pt', '1pt'),...
        FontSize('8pt'),...
        FontFamily('Arial')};    

    % clean up reactivity table - only show 3 sig digits for float data
    iso = F.out.Properties.VariableNames(...
        (find(strcmp(F.out.Properties.VariableNames, 'events')) + 1):end);      
    fltcols = [{'left', 'right'} iso];   
    O1 = varfun(@(x) num2str(x, ['%' sprintf('.%df', 3)]),...
        F.out(:, fltcols));
    O1.Properties.VariableNames = fltcols;    
    ha = {'name', 'size'};
    hb = [{'peak', 'left', 'right', 'events'}, iso];
    if usesplit == 0
        cols = {'size', 'peak', 'name', 'events'};
        hdr = [ha hb];
    else
        cols = {'size', 'bin', 'peak', 'name', 'events'};
        hdr = [ha {'bin'} hb];
    end
    O2 = [O1 F.out(:, cols)];
    tdata = [hdr; table2cell(O2(:, hdr))];
        
    % make Table object for reactivity data
    T = Table(tdata);
    T.TableEntriesStyle = {...
        VAlign('middle'),...
        HAlign('right'),...
        InnerMargin('2pt', '2pt', '2pt', '2pt'),...
        FontSize('9pt'),...
        FontFamily('Arial')};    
    T.Style = [T.Style, {...
        Width('100%'),...
        Border('solid'),...
        ColSep('solid'),...
        RowSep('solid')}];    
    T.Children(1).Style = {...
        VAlign('middle'),...
        HAlign('center'),...
        Bold,...
        BackgroundColor('lightblue')};

    % create report object
    R = Report([outdir fcs_base '_REPORT'], 'pdf');

    % overall organization
    add(R, TableOfContents);
    CH = struct;    
    tabs = {...
        'Table 1: Analysis Parameters and Count Summary';...
        'Table 2: Reactivity Data'};
    figs = {
        'SSC-H vs FSC-H to Identify Bead Populations';...
        'FSC-H vs FSC-A Singlet Gate';...
        'Split Gates';...
        'Peaks for Each Bead Population';...
        'Reactivities per Isotype'};
    if usesplit == 0, figs = figs([1 2 4 5]); end
    figs = cellstr("Figure " + num2str((1:length(figs))') + ": " + figs);
    tt = [tabs; figs];    
    for i = 1:(6 + usesplit)
        CH(i).obj = Chapter('Title', tt{i}, 'Numbered', false);
    end
    
    % profile settings
    add(CH(1).obj, {' '; ' '});
    add(CH(1).obj, M);
    add(CH(1).obj, {' '; ' '});
    add(CH(1).obj, PR);

    % reactivity table
    add(CH(2).obj, {' '; ' '});
    add(CH(2).obj, T);

    % figures    
    ifig = [1 2 3 4 5]; if usesplit == 0, ifig = [1 2 4 5]; end    
    for k = 1:length(ifig)
        img = Image(getSnapshotImage(Figure(F.figures(ifig(k))), R));
        img.Style = {ScaleToFit(true)};
        add(CH(k + 2).obj, {' '; ' '});
        add(CH(k + 2).obj, img);
    end
        
    % add chapters to report
    for i = 1:length(CH), add(R, CH(i).obj); end

    % close report
    close(R)
        
