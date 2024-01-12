%
% For a linear classifier (SVM, LDA), define the hyperplane and the signal axis
%
% Beta/Bias: parameters as derived from fitclinear or fitcdiscr
% P: Dimensionality of the data
%
%urut/Aug'23
function [signalVector_derived, signalVector_fromModel_norm, fHyperplane, fBoundLine] = LinearClassifier_getSignalAxis(Betas, Bias, P )

%Beta = mdl.Beta;
%b    = mdl.Bias;

Beta = Betas;
b    = Bias;

%get arbitrary points that fall on the hyperplane. x1 is free parameters, x2...xn are set to arbitrary values
fHyperplane = @(X) -(X(2:end)*Beta(2:end) + b)/Beta(1);
fBoundLine = @(X) -(X*Beta(1) + b)/Beta(2);

%get P-1 random pairs of points on the hyperplane
H=[];
for k=1:P-1    
    %point 1 in pair
    hypPoint1 = randn(1,P);
    hypPoint1(1) = fHyperplane(hypPoint1);

    %point 2 in pair
    hypPoint2 = randn(1,P);
    hypPoint2(1) = fHyperplane(hypPoint2);

    H(:,k) = hypPoint2-hypPoint1;
end

%== derive signal vector from hyperplane using QR decomposition
% follows: https://math.stackexchange.com/questions/3810098/how-do-you-get-an-orthogonal-vector-in-arbitrary-dimensions
[Q,R] = qr(H);
l = size(H,2);   % l+1, l+2 ... are orthogonal
signalVector_derived = Q(:,l+1);

%check - this should all be orthogonal to all points and norm=1
angels=[];
for k=1:l
    angels(k) = dot(signalVector_derived, H(:,k)/norm(H(:,k)) );
end

if ~isempty(find(angels>0.001))
    warning('signalVector is not orthogonal to hyperplane, check manual');
end

%Define signal vector directly from weights to see that it is parallel
signalVector_fromModel = Beta;
signalVector_fromModel_norm = signalVector_fromModel/norm(signalVector_fromModel);

%Should should be =1,-1
angl_sigVectors = dot( signalVector_fromModel_norm, signalVector_derived);

if abs(angl_sigVectors)>1.0001 || abs(angl_sigVectors)<0.0009
    warning('signalVector is not parallel to model parameters, check manual');
end


