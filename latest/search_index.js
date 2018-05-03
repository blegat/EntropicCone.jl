var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SumOfSquares-–-Entropic-Cone-1",
    "page": "Index",
    "title": "SumOfSquares –- Entropic Cone",
    "category": "section",
    "text": "EntropicCone utilities for constructing and using outer and inner approximations of the Entropic Cone."
},

{
    "location": "index.html#Contents-1",
    "page": "Index",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"intro.md\"]\nDepth = 2"
},

{
    "location": "intro.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "intro.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": ""
},

{
    "location": "intro.html#The-definition-of-entropy-1",
    "page": "Introduction",
    "title": "The definition of entropy",
    "category": "section",
    "text": "In 1948, Shannon published \"A Mathematical Theory of Communication\" [Sha48]. In this paper, Shannon introduces the entropy of a random variable. Suppose we have a random variable X of alphabet mathcalX, he defines the entropy of X asH_b(X) = sum_x in mathcalX PrX = x log_b frac1PrX = xwhere the basis b is positive. If b is 2 (resp. e), the unit is the bits (resp. nats). Note that H_b(X) = H_a(X) log_b(a) so the entropies using different basis are equivalent up to a positive constant factor.The entropy of several random variables in simply the entropy of their cartesian product:beginalign*\n  H_b(X_1 ldots X_n)\n   = sum_(x_1ldotsx_n) in mathcalX_1 times cdots times mathcalX_n Pr(X_1 ldots X_n) = (x_1ldotsx_n) log_b frac1Pr(X_1 ldots X_n) = (x_1ldotsx_n)\nendalign*By convention, we say that the entropy of an empty set of random variables is 0.Given a n random variables, we can compute the entropy of any of the 2^n subset of those n variables. The entropic vector of a set of n random variables is a vector h, indexed by the subsets of n = 1 ldots n, such that h_S = H_b( X_i mid i in S)."
},

{
    "location": "intro.html#The-entropic-cone-1",
    "page": "Introduction",
    "title": "The entropic cone",
    "category": "section",
    "text": "The entropic cone of n variables is the set of vectors of mathbbR^2^n-1 that are entropic:mathcalH_n =  h in mathbbR^2^n-1 mid exists X_1 ldots X_n forall emptyset neq S subseteq n h_S = H_b( X_i mid i in S) We do not include the dimension corresponding to the entropy of the empty set as it is zero to make the cone mathcalH_n solid, i.e. full-dimensional.[Sha48] Claude Elwood Shannon. A mathematical theory of communication. Bell System Technical Journal, 27:379–423 and 623–656, July and October 1948."
},

]}
