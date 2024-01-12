function [pval,Factual,Fdist] = randanova1(m,group,num)

% function [pval,val,dist] = randANOVA1(m,group,num,mode)
%
% Estimate p-value for one-way ANOVA using permutation.
% For further information on using permutation to calculate signficance
% tests see Anderson (doi:10.1139/cjfas-58-3-626).
%
% <m> is a vector to be permuted
% <group> is grouping variable
% <num> (optional) is number of randomizations.  default: 1000.
% 
% return p-values in <pval>, actual F value in <Factual>, and 
% distribution of randomly obtained F values in <Fdist>.
%
% Can be run in parallel mode(type: matlabpool open)
%
% e.g.
%   y = [52.7 57.5 45.9 44.5 53.0 57.0 45.9 44.0]';
%   g1 = [1 2 1 2 1 2 1 2];
%   g2 = {'hi';'hi';'lo';'lo';'hi';'hi';'lo';'lo'};
%
%   [pval,Factual,~]=randanova1(y,g1)
%
%   pval =
%
%       0.6520
%
%
%   Factual =
%
%        0.1043
%
%   [pval,Factual,~]=randanova1(y,g2)
%
%   pval =
%
%       0.0290
%
%
%   Factual =
%
%        53.3575
%
%
%
% $	Author: David Stern	$   $   Date :2013/07/30   $
% $ Janelia Farm Research Campus, HHMI, Ashburn, VA
%


% inputs
if ~exist('num','var') || isempty(num)
  num = 1000;
end
fun = (@(x) anova1(x,group,'off'));

% calc actual
[~,Table,~] = feval(fun,m);
Fcolumn = find(strcmp('F',Table(1,:)));
Factual = cell2mat(Table(2,Fcolumn));

%make randomized matrix of num columns of permuted data
Fdist = zeros(num,1);
parfor p=1:num
    permm = datasample(m,numel(m),'Replace',false);
    [~,Table,~] = feval(fun,permm);
    Fdist(p) = cell2mat(Table(2,Fcolumn));
end

pval = sum(ge(Fdist,Factual)) / num;