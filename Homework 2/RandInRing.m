function [A, B] = RandInRing(r, w, d, n)

% Build A region
Xrand = 2*(r + w)*rand(n,1) - (r + w);
Yrand = (r + w)*rand(n,1);
randMat = [Xrand, Yrand];
ValidCoords = 0;
while ValidCoords < n
    ii = 1;
    while ii <= length(randMat)
        if (norm(randMat(ii,:)) > r + w) || (norm(randMat(ii,:)) < r)
            randMat(ii,:) = [];
            ii = ii - 1; 
        end
        ii = ii + 1;
    end
    ValidCoords = length(randMat);
    Xrand = 2*(r + w)*rand(n - ValidCoords,1) - (r + w);
    Yrand = (r + w)*rand(n - ValidCoords,1);
    randMat = [randMat; [Xrand, Yrand]]; %#ok
end
A = randMat;

% Build B region
Xrand = 2*(r + w)*rand(n,1) - w/2;
Yrand = -(d + (r + w)*rand(n,1));
randMat = [Xrand, Yrand];
ValidCoords = 0;
while ValidCoords < n
    ii = 1;
    while ii <= length(randMat)
        if (norm(randMat(ii,:) - [r + w/2, -d]) > r + w) || ...
                (norm(randMat(ii,:) - [r + w/2, -d]) < r)
            randMat(ii,:) = [];
            ii = ii - 1;
        end
        ii = ii + 1;
    end
    ValidCoords = length(randMat);
    Xrand = 2*(r + w)*rand(n - ValidCoords,1) - w/2;
    Yrand = -(d + (r + w)*rand(n - ValidCoords,1));
    randMat = [randMat; [Xrand, Yrand]]; %#ok
end
B = randMat;
end

