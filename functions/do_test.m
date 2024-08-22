function [d] = do_test(a,varargin)
defaultbp = false;
defaultbpwidth = [150];

p = inputParser;
validScalarPosNum = @(x) isnumeric(x);
validbool = @(x) islogical(x);
addRequired(p,'a',validScalarPosNum);
addOptional(p, 'bp', defaultbp, validbool);
addOptional(p,'bpwidth',defaultbpwidth,validScalarPosNum);
parse(p,a,varargin{:});

bpwidth = p.Results.bpwidth;
if p.Results.bp == true
    d = a + bpwidth
elseif defaultbp == false
    d = a;
end 

end 
