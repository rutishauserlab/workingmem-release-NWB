% modeNoiseAxis: 1 use PCA, 2 use covariance matrix of entire data, 3 use average noise covariance matrix
%   Node for modeNoiseAxis=2: this will also include signal correlations. This is useful when comparing 
%   against a within-group scramble, in which signal corrs will be preserved.
%
%urut/aug23
function [errRate, signalVector_toUse, noiseVector1, noiseVector2, meanAngle, projValues_Cl1, projValues_Cl2, fHyperplane, fBoundLine, PP] ...
    = popGeometry_defineAxes(modeClassifier, data, labels, Pdims, NperGroup, modeNoiseAxis)

switch(modeClassifier)
    case 1      
        obj_fisher1 = fitcdiscr(data, labels);
        cvmodel = crossval(obj_fisher1);
        errRate = kfoldLoss(cvmodel)

        K = obj_fisher1.Coeffs(2).Const;
        L = obj_fisher1.Coeffs(2).Linear;

        [signalVector_derived, signalVector_fromModel_norm, fHyperplane] = LinearClassifier_getSignalAxis(L, K, Pdims);
        signalVector_toUse = signalVector_derived;
        
    case 2
        obj_svm1 = fitclinear( data, labels,'KFold',5);
        errRate = kfoldLoss(obj_svm1);    % to see problem is solvable
        mdl = obj_svm1.Trained{1};  % pick one of the trained models to plot

        [signalVector_derived, signalVector_fromModel_norm, fHyperplane, fBoundLine] = LinearClassifier_getSignalAxis(mdl.Beta, mdl.Bias, Pdims);
        signalVector_toUse = signalVector_derived;
end

indsGrp1 = find(labels==0);
indsGrp2 = find(labels==1);

switch(modeNoiseAxis)
    case 1         %Noise vectors using PCA, separately for each class


        % Define noise vectors
        [coeff_Cl1] = pca( data( indsGrp1,:) );
        noiseVector1 = coeff_Cl1(:,1);

        [coeff_Cl2] = pca( data( indsGrp2,:) );
        noiseVector2 = coeff_Cl2(:,1);

        %Angle. https://stackoverflow.com/questions/18330451/angle-between-two-vectors-matlab
        angle1 = subspace(signalVector_toUse, noiseVector1 );
        angle2 = subspace(signalVector_toUse, noiseVector2 );

        meanAngle = mean([angle1 angle2]);
    case 2         %Noise vectors using covariance matrix of the entire data (pooled)

        %use covariance matrix to define noise axis
        CorrMatrix = cov(data);

        % compute noise axis based on covariance matrix.  Following Averbeck & Lee 2006, J neurophysiol
        [V,D] = eig( CorrMatrix );
        indMaxEig = find( diag(D) == max(diag(D)));
        indMaxEig = indMaxEig(1);

        noiseVector1 = V(:,indMaxEig);  % eigenvector associated with max eigenvalue
        noiseVector2 = noiseVector1;  %they are the same

        meanAngle = subspace(signalVector_toUse, noiseVector1);

    case 3  %Noise vectors using covariance matrix of each class

        CorrMatrix1 = cov( data(indsGrp1,:) );
        CorrMatrix2 = cov( data(indsGrp2,:) );

        CorrMatrix_Mean = (CorrMatrix1+CorrMatrix2)./2;

        [V,D] = eig( CorrMatrix_Mean );
        indMaxEig = find( diag(D) == max(diag(D)));
        indMaxEig = indMaxEig(1);

        noiseVector1 = V(:,indMaxEig);  % eigenvector associated with max eigenvalue
        noiseVector2 = noiseVector1;  %they are the same

        meanAngle = subspace(signalVector_toUse, noiseVector1);


end

%==
%Projected Precision (PP). Following Nogueira et al 2020 J Neurosic and Roussy et al 2021 Molecular Psychiatry
%PP is a weighted sum of all the noise axis angles.
%
PP=0;
if modeNoiseAxis==2 || modeNoiseAxis==3
    % A is the covariance matrix
    if modeNoiseAxis==2
        A = CorrMatrix;
    else
        A = CorrMatrix_Mean;
    end

    if ~issymmetric(A)
        warning('covariance not symmetric');
    end

    %== use SVD to derive orthonormal basis.  
    %A = U*S*V'
    %U==V because A is symmetric
    %A = U*S*U'
    %U(:,k) are the eigenvectors, which are orthonormal
    %diag(S) are the singular values, which are refererd to as sigma_k^2 in below papers
    [U,S,V] = svd(A);
    eigVectors = U;
    singularValues = diag(S);

    %loop over all noise vectors
    PP=0;
    ang_k=[];
    for k=1:length(singularValues)
        %Eq 4 in Nogueira 2020;   See Eqs 71-76 in supplement of Moreno-Bote et al 2017 Nat Neurosci for how this is derived.
        % The assumption is that the eigVectors(:,k) are orthonormal
        ang_k(k) = subspace(signalVector_toUse, eigVectors(:,k));   
        PP = PP + cos( ang_k(k) )^2/singularValues(k);
    end

    anglSum = sum( cos(ang_k).^2 );   % sanity check - should be equal to 1. See Eq 76 in supplement of Moreno-Bote et al 2017 Nat Neurosci
    if anglSum <0.99 | anglSum>1.01
        warning(['sum of square angles not equal=1? is=' num2str(anglSum)]);
    end

    PP = sqrt(PP);
else
    %cant be computed with other methods
    PP=0;
end


% project the data onto the signal axis
projValues_Cl1 = getNoiseVar_alongSignalAxis ( data(1:NperGroup,:), signalVector_derived);
projValues_Cl2 = getNoiseVar_alongSignalAxis ( data(NperGroup+1:end,:), signalVector_derived);

