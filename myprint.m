function[] = myprint(varargin)

writeoption = varargin{1};
file        = varargin{2};
text        = varargin{3};
variables   = varargin{4:end};

if strcmp(writeoption,'Y')   
    fprintf(file,text,variables);
end

end