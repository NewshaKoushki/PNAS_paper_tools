function m = dlmread(filename,dlm,r,c,rng)
%DLMREAD Read ASCII delimited file.
%   M = DLMREAD(FILENAME,DLM) reads numeric data from the ASCII delimited
%   file FILENAME using the delimiter DLM.  The result is returned in M.
%   Use '\t' to specify a tab.
%
%   M = DLMREAD(FILENAME,DLM,R,C) reads data from the DLM-delimited
%   file FILENAME.  R and C specify the row R and column C
%   where the upper-left corner of the data lies in the file.  R and C
%   are zero-based so that R=0 and C=0 specifies the first value in the file.
%
%   M = DLMREAD(FILENAME,DLM,RNG) reads the range specified
%   by RNG = [R1 C1 R2 C2] where (R1,C1) is the upper-left corner of
%   the data to be read and (R2,C2) is the lower-right corner.  RNG
%   can also be specified using spreadsheet notation as in RNG = 'A1..B7'.
%
%   DLMREAD fills empty delimited fields with zero.  Data files where
%   the lines end with a non-space delimiter will produce a result with
%   an extra last column filled with zeros.
%
%   See also DLMWRITE, TEXTREAD, WK1READ, WK1WRITE.

% Obsolete syntax:
%   M = DLMREAD(FILENAME,DLM,R,C,RNG) reads only the range specified
%   by RNG = [R1 C1 R2 C2] where (R1,C1) is the upper-left corner of
%   the data to be read and (R2,C2) is the lower-right corner.  RNG
%   can also be specified using spreadsheet notation as in RNG = 'A1..B7'.
%   A warning will be generated if R,C or both don't match the upper
%   left corner of the RNG.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.26 $  $Date: 1998/12/10 22:07:12 $

offset = 0; % True if data written must be offset

% Parse inputs
if nargin<1, error('Not enough input arguments.'); end
if ~isstr(filename), error('Filename must be a string.'); end

if nargin==1 % delimiter defaults to Comma for CSV
  dlm = ',';
else
  dlm = sprintf(dlm); % Interpret \t (if necessary)
end
if length(dlm) ~= 1,
  error('DLM must be a single character.');
end
 
if nargin<=2, % dlmread(file) or dlmread(file,dim)
  r = 0;
  c = 0;
  nrows = -1; % Read all rows
elseif nargin==3, % dlmread(file,dlm,rng)
  if length(r)==1, % Catch obsolete syntax dlmread(file,dlm,r)
    warning('Obsolete syntax. C must be specified with R.');
    m = dlmread(filename,dlm,r,0);
    return
  end
  rng = r;
  if isstr(rng)
    rng = str2rng(rng);
  end
  r = rng(1);
  c = rng(2);
  nrows = rng(3) - rng(1) + 1;
elseif nargin==4, % dlmread(file,dlm,r,c)
  nrows = -1; % Read all rows
elseif nargin==5, % obsolete syntax dlmread(file,dlm,r,c,rng)
  if isstr(rng)
    rng = str2rng(rng);
  end
  rold = r; cold = c; offset = 1;
  if r > rng(3) | c > rng(4), m = []; return, end
  if r ~= rng(1) | c ~= rng(2)
    warning('R and C should match RNG(1:2).  Use DLMREAD(FILE,DLM,RNG) instead.')
  end
  % For compatibility
  r = max(rng(1),r);
  c = max(rng(2),c);
  nrows = rng(3) - r + 1;
end

% Determine the number of columns in the file 
% by looking at the firstline
fid = fopen(filename,'r');
if (fid == -1)
  error('File not found or permission denied.');
end
firstline = fgetl(fid);
fclose(fid);
if isempty(firstline), error('File is empty.'); end

if dlm == ' '
  count = length(findrun(firstline == dlm));
  if firstline(1)==' ', count = count-1; end
else 
  count = sum(firstline == dlm);
end
ncols = count + 1;

% Determine format for textread and check for common error conditions
if nargin==3 | nargin==5
  n = rng(4) - c + 1;
  if rng(4) >= ncols, error('Range outside data.'); end
  if n < 1, error('Column C outside data.'); end
  format = [repmat('%*s',1,c) repmat('%f',1,n)  repmat('%*s',1,ncols-c-n)];
else
  n = ncols - c;
  if n < 1, error('Column C outside data.'); end
  format = [repmat('%*s',1,c) repmat('%f',1,n)];
end

% Try reading the file using textread since it is faster.  If we
% get an error make note of it and try old_dlmread for compatibility.
try
  errflg = 0;
  lasterr('')
  % Read the file using textread
  if dlm == ' '
    [cols{1:n}] = textread(filename,format,nrows,'headerlines',r);
  else
    dlm = sprintf([dlm '\n']);
    whitespace = setdiff(sprintf(' \b\r\t'),dlm);
    [cols{1:n}] = textread(filename,format,nrows,'delimiter',dlm, ...
                     'whitespace',whitespace,'headerlines',r);
  end

  % Special case -- non-space delimiters at the end of a line and
  % last element wasn't read.
  k = length(deblank(firstline));
  if firstline(k) == dlm(1) & length(cols{end}) ~= length(cols{1}) & ...
     all(cols{end}==0)
    cols{end} = [cols{end};0];
  end
catch
  errflg = 1;
  errmsg = lasterr;
  cols = {};
end

lengths = cellfun('length',cols);
if ~errflg & all(lengths == lengths(1))
  m = [cols{:}];
else
  % Try reading this file with the old dlmread -- for compatibility
  offset = 0; % Disable offset effect
  try
    switch nargin
    case {0,1}
      m = old_dlmread(filename);
    case 2
      m = old_dlmread(filename,dlm(1));
    case 3
      m = old_dlmread(filename,dlm(1),r,c,rng);
    case 4
      m = old_dlmread(filename,dlm(1),r,c);
    case 5
      m = old_dlmread(filename,dlm(1),r,c,rng);
    end
  catch
    error(errmsg)
  end
end

% For backwards compatibility
if nargin==5 & offset
  mm((r+1:rng(3)+1),(c+1:rng(4)+1)) = m;
  m = mm(rold+1:end,cold+1:end);
end

%------------------------------
function [b,e] = findrun(v)
%FINDRUN Find runs of non-zero elements in a vector.
%   [B,E] = FINDFUN(V) returns the indices of the beginning B
%   and end E of the runs of non-zero elements in the vector V.
%
%   Examples:
%     [b,e] = findrun([0 1 1 0 1 1 0 3 1]) returns 
%         b = [2 5 8] and e = [3 6 9]

d = diff([0 v(:)'~=0 0]);
b = find(d==1);
e = find(d==-1)-1;

%-----------------------------
function m=old_dlmread(filename, dlm, r, c, rng)
%DLMREAD Read ASCII delimited file.
%   This is the 5.2 version of dlmread and is provided as this
%   subfunction for compatibility.  This will only be called
%   if the current version of dlmread fails.
%
%   M = DLMREAD(FILENAME,DLM) reads the data from the ASCII delimited
%   file FILENAME using the delimiter DLM.  The result is returned in M.
%   Use '\t' to specify a tab. The file can only contain numeric values.
%
%   M = DLMREAD(FILENAME,DLM,R,C) reads data from the DLM-delimited
%   file format FILENAME.  R and C specify the row R and column C
%   where the upper-left corner of the data lies in the file.  R and C
%   are zero-based so that R=0 and C=0 specifies the first value in the file.
%
%   M = DLMREAD(FILENAME,DLM,R,C,RNG) reads only the range specified
%   by RNG = [R1 C1 R2 C2] where (R1,C1) is the upper-left corner of
%   the data to be read and (R2,C2) is the lower-right corner.  RNG
%   can also be specified using spreadsheet notation as in RNG = 'A1..B7'.

% test for proper filename
if nargin<1, error('Not enough input arguments.'); end
if ~isstr(filename), error('Filename must be a string.'); end

NEWLINE = sprintf('\n');

% check/set row,col offsets
if nargin<3, r = 0; end
if nargin<4, c = 0; end

% delimiter defaults to Comma for CSV
if nargin<2, dlm = ','; end
dlm = sprintf(dlm); % Handles special characters.

% get the upper-left and bottom-right cells of the range to read into MATLAB
if nargin==5
    if ~isstr(rng)
        ulc = rng(1:2);
        brc = rng(3:4);
    else
        x = str2rng(rng);
        ulc = x(1:2);
        brc = x(3:4);
    end
    all = 0;
else
    all = 1;
    rng = [ 0 0 ];
    ulc = [0 0];
    brc = [0 0];
end

ulc = fliplr(ulc+1);
brc = fliplr(brc+1);

% open the file 
fid = fopen(filename,'r');
if fid == (-1)
    error(['dlmread: Could not open file filename ']);
end

% Read delimited format 
eol = NEWLINE;    % End Of Line char
loc = [1 1];      % starting location of return matrix
line = fgets(fid); % get the 1st line, if any...

% read till eof
while ~isequal(line,-1)
    i = 1;
    j = 1;
    while i <= length(line)
        %
        % read chars from line, parsing delimiters & numbers
        %
        num = [];
        j = 1;

        while (i <= length(line)) & (line(i) ~= dlm) & (line(i) ~= eol)
            % build number string from characters on the line
            num(j) = line(i);
            i = i + 1;    % overall line index
            j = j + 1;    % number string index
        end

        % found a delimiter or <eol>
        if(all | ((loc >= ulc) & (loc <= brc)))
            if ~isempty(num)
                v=str2num(setstr(num));
                if ~isempty(v)
                    m(loc(2)+r, loc(1)+c) = v;
                else % try a char
                    m(loc(2)+r, loc(1)+c) = setstr(num);
                end
            else     
                % no number found between delimiters
                m(loc(2)+r, loc(1)+c) = 0;
            end
        end

        if (i <= length(line)) & (line(i) == dlm)
            % delimiter, set location to next row and get next line
            loc(1) = loc(1) + 1;
            i = i + 1;
        else
            if (i > length(line)) | (line(i) == eol)
                % eol, set location to next row and get next line
                loc(2) = loc(2) + 1;
                loc(1) = 1;
                i = i + 1;
            end
        end
    end
    % get next line of file
    line = fgets(fid); 
end
% close file
fclose(fid);

% Return only the valid part
m = m(2*r+1:end,2*c+1:end);
