function opt = parseOptArgs (defaults, varargin)
%
% PARSEOPTARGS - Optional arguments parsing
%   
% SYNTAX
%
%   OPT = PARSEOPTARGS( DEFAULTS, 'Name', Value, ... )
%
% INPUT
%
%   DEFAULTS    Struct with default parameters          [struct]
%               (field names and values comprise the
%               recognized name-value pairs)
%   <Name-value pairs for non-default parameters>       [varargin]
%(1:7)
% OUTPUT
%
%   OPT         Input and default parameters            [struct]
%
% DESCRIPTION
%
%   OPT = PARSEOPTARGS(DEFAULTS,'Name',Value,...) is a wrapper around
%   MATLAB's inputParser, to facilitate optional argument parsing to
%   functions.
%
%   If any value is empty (e.g., []), then the default value is used.
%
% PARAMETERS
%
%   The following parameters are explicitly set for the underlying
%   inputParser.  Use a customized argument parser if these are not
%   appropriate.
%
%       CaseSensitive   = false
%       KeepUnmatched   = false
%       PartialMatching = true
%       StructExpand    = true
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      inputParser
%
    
    
    %% INITIALIZATION
    
    ip = inputParser;
    
    ip.CaseSensitive   = false;
    ip.KeepUnmatched   = true;
    ip.PartialMatching = true;
    ip.StructExpand    = true;
    
    argNames = fieldnames( defaults );
    
    
    %% PARAMETERS
    
    for i = 1 : length(argNames)
        addParameter( ip, argNames{i}, defaults.(argNames{i}) );
    end
    
    
    %% PARSE AND RETURN
    
    parse( ip, varargin{:} );
    
    opt = ip.Results;
    
    
    %% SET EMPTY VALUES TO DEFAULTS
    
    for i = 1 : length(argNames)
        if isempty( opt.(argNames{i}) )
            opt.(argNames{i}) = defaults.(argNames{i});
        end
    end
    
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% VERSION
%
%   0.1 - April 21, 2017
%
% CHANGELOG
%
%   0.1 (Apr 21, 2017) - Alexandros
%       * initial implementation
%
% ------------------------------------------------------------
