% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/butte.R
\name{.lpbounds}
\alias{.lpbounds}
\title{Use linear programming to solve the bounds of TK or T0}
\usage{
.lpbounds(q, possible_histories, scost = 100, p0 = FALSE)
}
\arguments{
\item{q}{q estimated from data}

\item{possible_histories}{matrices of possible SCNA-SSNV histories, see also function "cnmutHistory"}

\item{scost}{the cost for slack variables (default 100)}

\item{p0}{logical, if TRUE, the upper bounds for T0 will be estimated (instead of TK)}
}
\value{
the lower and upper bounds of the time duration for the last stage
}
\description{
This function estimate the bounds of TK, which is the last time interval in the tumor cell envolution
history. Based on different possible history given a copy number ratio, the function makes use of 
linear programming to minimize and maximize TK. The objective function of the optimization if f(x) = ta, which
can be written as [0 0 0 0 0 ... 0 1]T*[t1 t2 t3 ... ta]. The first part of constraints are given by (A-qs')t = 0, 
where s' refers to the transpose vector calculated by rowSum(A). The second part of constraints is the 
convexity of time vector t. Each element of t refers to a relative fraction of time. 
We then combine these two constraints into a single linear system formation. he first constraint directly 
follows from At/c = q, wherer c is a normalizing constant given by the product s'*t.
}
\details{
Then the feasibility of the solution region of the linear programming problem will only be influenced by q, which
is optimized previously from data. The uncertainty of q can make the solution space infeasible. So we add some 
slack variables to elasticize the linear programming problem. 
For details, please check: http://web.mit.edu/lpsolve/doc/Infeasible.htm
This elasticizing method will find the approximate bounds of TK close to the constraints.
"scost" is the argument adjusting the penalty of the additional slack variables.
}
