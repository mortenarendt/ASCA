function Xi = getContinousHat(y,p)

%%% a bunch of function forms

id = ceil(rand(p,1)*5); 

for i=1:p
    idd= id(i);
    switch idd
        case 1
            % sinusoides
            x = getSinus(y);
        case 2
            % inverse
            e = 1;
            x = randn(1)./ (y - min(y) + e);
        case 3
            % exponential
            yn = y / range(y);
            x = exp(randn(1)*yn + randn(1));
            
        case 4
            % linear
            x = randn(1)*y + randn(1);
        case 5
            % polinomial
            x = getPoly(y,4);
    end
    x = (x - mean(x)) / std(x); 
    Xi(:,i) = x; 
end


function x = getSinus(y)
r = rand(1);
k = 1/range(y)*pi*(r+0.5);
sh = range(y) + rand(1);
x = sin(y * k + sh);

function x = getPoly(y,mxorder)

ord = ceil( rand(1)*mxorder);
B = [ones(length(y),1)];
for i=1:ord
    B(:,i+1) = y.^i;
end
x = B*randn(ord+1,1);
