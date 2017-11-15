function[] = myprint(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to fix printing of messages and output to external files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writeoption = varargin{1};
file        = varargin{2};
text        = varargin{3};

if strcmp(writeoption,'Y')  
    if nargin == 3 % Equivalent to disp command, no variable input
        fprintf(file,text);
    elseif nargin > 3 % Printing of variables
        variables   = varargin{4:end};
        fprintf(file,text,variables);
    end
end

end