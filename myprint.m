function[] = myprint(varargin)

writeoption = varargin{1};
file        = varargin{2};
variables   = varargin{3:end};

if strcmp(writeoption,'Y')   
    fprintf(file,text,variables);
end

end