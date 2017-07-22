function[] = myprint(varargin)

writeoption = varargin{1};
file        = varargin{2};
text        = varargin{3};

if strcmp(writeoption,'Y')  
    if nargin == 4
        fprintf(file,text);
    else
        variables   = varargin{4:end};
        fprintf(file,text,variables);
    end
end

end