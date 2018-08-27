function generateData()
   y1 = [
    0.27  0.33 ;
    0.35  0.43 ;
    0.42  0.5 ;
    0.35  0.61;
    0.29  0.74;
    0.49  0.6 ;
    0.49  0.43;
    0.58  0.34;
    0.59  0.75;
   ];
   y2 = [
   0.15   0.45;
   0.155  0.55;
   0.39   0.5 ;
   0.52   0.32;
   0.5    0.5 ;
   0.51   0.73;
   0.62   0.61;
   0.62   0.42;
   0.76   0.5 ;
   ];
   Y1 = generateMeasures(y1);
   Y2 = generateMeasures(y2);
   Y1, Y2
end


function generateMeasures(X)
    nMeasures = 10;
    a = Array{Float64}[];
    for i in 1:nMeasures
        push!(a,generateNoisyEstimate(X));
    end
    a
end

function generateNoisyEstimate(X)
    rmvec = []
    for i in 1:size(X)[1]
        if rand() < 1/(4*size(X)[1])
            #deleteat!(X)
            rmvec = [rmvec; i];
        end
    end
    X = X[filter(x->!(x in rmvec),1:size(X)[1]),:];
    for i in 1:size(X)[1]
        if rand() < 1/2
            X[i,:] = randomizePosition(X[i,:]);
        end
    end
    X
end

function randomizePosition(x)
    x[1] += randn()/100;
    x[2] += randn()/100;
    x
end

