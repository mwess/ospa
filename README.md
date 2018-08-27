# ospa

[![Build Status](https://travis-ci.org/mwess/ospa.jl.svg?branch=master)](https://travis-ci.org/mwess/ospa.jl)

[![Coverage Status](https://coveralls.io/repos/mwess/ospa.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/mwess/ospa.jl?branch=master)

[![codecov.io](http://codecov.io/github/mwess/ospa.jl/coverage.svg?branch=master)](http://codecov.io/github/mwess/ospa.jl?branch=master)


Implementation of the ospa barycenter[[1]](https://ieeexplore.ieee.org/document/7266717/).

Additionally two initialization procedures are added: Mean-initialization and iterative-initialization. 

Usage: 

`using ospa`


Examples:

Example Data can be loaded:

`X,Y = ospa.generateData()`

Using different initialization methods:

Mean initialization method:

`ospa_barycenter(X, 0.05, missing, :binned)`

Iterative initialiation methods:

`ospa_barycenter(X, 0.05, missing, :iterative)`

Standard initialization (initialize cluster a random measurement):

`ospa_barycenter(X, 0.05, missing, :standard)`
