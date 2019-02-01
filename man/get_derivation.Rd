\name{get_derivation}
\alias{get_derivation}
\title{Identify a causal effect from arbitrary experiments and observations}
\description{Identify a causal \code{query} from available \code{data} in a causal model described by a semi-Markovian \code{graph}. Special mechanisms related to transportability of causal effects, recoverability from selection bias and identifiability under missing data can also be included in the model. }

\usage{get_derivation(data, query, graph, transportability = NULL, selection_bias = NULL,
  missing_data = NULL, control = list())}

\arguments{
  \item{data}{a character string describing the available distributions in the package syntax. See \sQuote{Details}}
  \item{query}{a character string describing the target distribution in the package syntax. See \sQuote{Details}}
  \item{graph}{a character string describing the graph in the package syntax. See \sQuote{Details}}
  \item{transportability}{a character string describing the transportability nodes of the model in the package syntax. See \sQuote{Details}}
  \item{selection_bias}{a character string describing the selection bias nodes of the model in the package syntax. See \sQuote{Details}}
  \item{missing_data}{a character string describing the missing data mechanisms of the model in the package syntax. See \sQuote{Details}}
  \item{control}{a \link{list} of control parameters. See \sQuote{Details}. }
}

\details{

\code{data} is used to list the available input distributions of the form

\ifelse{html}{\out{<p> P(A<sub>i</sub>|do(B<sub>i</sub>),C<sub>i</sub>) </p>}}{\deqn{P(A_i|do(B_i),C_i).}}

Individual variables within sets should be separated by a comma. For example, three input distributions

\ifelse{html}{\out{<p> P(Z|do(X)), P(W,Y|do(Z,X)), P(W,Y,X|Z) </p>}}{\deqn{P(Z|do(X)), P(W,Y|do(Z,X)), P(W,Y,X|Z),}}

should be given as follows: 
\preformatted{
> data <- "
+  P(Z|do(X))
+  P(W,Y|do(Z,X))
+  P(W,Y,X|Z)
+"
}

The use of multiple do-operators is not permitted. Furthermore, when both conditioning variables and a do-operator are present, every conditioning variable must either precede the do-operator or follow it.

\code{query} is the target distribution of the search. It has the same syntax as \code{data}, but only a single distribution should be given.

\code{graph} is a description of a directed acyclic graph where directed edges are denoted by \code{->} and bidirected arcs corresponding to unobserved confounders are denoted by \code{--}. As an example a graph with two directed edges and one bidirected arc is constructed as follows: 
\preformatted{
> graph <- "
+  X -> Z
+  Z -> Y
+  X -- Y
+"
}

\code{transportability} enumerates the nodes that should be understood as transportability nodes responsible for discrepancies between domains. Individual variables should be separated by a comma. See e.g., Bareinboim and Pearl (2014) for details on transportability.

\code{selection_bias} enumerates the nodes that should be understood as selection bias nodes responsible for bias in the input data sets. Individual variables should be separated by a comma. See e.g., Bareinboim and Pearl (2014) for details on selection bias recoverability.

\code{missing_data} enumerates the missingness mechanisms of the model. The syntax for a single mechanism is \code{M_X : X} where \ifelse{html}{\out{M<sub>X</sub>}}{\eqn{M_X}} is the mechanism for \eqn{X}. Individual mechanisms should be separated by a comma. Note that both \ifelse{html}{\out{M<sub>X</sub>}}{\eqn{M_X}} and \eqn{X} must be present in the graph if the corresponding mechanism is given as input. Proxy variables should not be included in the graph, since they are automatically generated based on \code{missing_data}. See e.g., Mohan, Pearl and Tian (2013) for details on missing data as a causal inference problem.

The \code{control} argument is a list that can supply any of the following components:
\describe{
  \item{\code{benchmark}}{A logical value. If \code{TRUE}, the search time is recorded and returned (in milliseconds). Defaults to \code{FALSE}. }
  \item{\code{draw_derivation}}{A logical value. If \code{TRUE}, a string representing the derivation steps as a dot graph is returned. Defaults to \code{FALSE}. }
  \item{\code{draw_all}}{A logical value. If \code{TRUE} and if \code{draw_derivation = TRUE}, the derivation will contain every step taken by the search. If \code{FALSE}, only steps that resulted in an identifiable target are returned. Defaults to \code{FALSE}. }
  \item{\code{formula}}{A logical value. If \code{TRUE}, a string representing the identifiable query is returned. If \code{FALSE}, only a logical value is returned that takes the value \code{TRUE} for an identifiable target and \code{FALSE} otherwise. Defaults to \code{TRUE}. }
  \item{\code{heuristic}}{A logical value. If \code{TRUE}, new distributions are expanded according to a search heuristic (see the vignette for details). Otherwise, distributions are expanded in the order in which they were identified. Defaults to \code{TRUE} unless missing data mechanisms are provided in \code{missing_data}. }
  \item{\code{improve}}{A logical value. If \code{TRUE}, various improvements are applied to the search to make it more efficient (see the vignette for details). Defaults to \code{TRUE}. }
  \item{\code{md_sym}}{A single character describing the symbol to use for active missing data mechanisms. Defaults to \code{"1"}. }
  \item{\code{replace}}{A logical value. If \code{TRUE}, the search will continue deriving new distributions after the target has been reached in order to possibly find a shorter search path. Defaults to \code{FALSE}. }
  \item{\code{rules}}{A numeric vector describing a subset of the search rules to be used. Must be a subset of \code{c(-1,-2,-3,1,2,3,4,5,-6,6,-7,7,8,9)} (see the vignette for details). }
  \item{\code{verbose}}{A logical value. If \code{TRUE}, diagnostic information is printed to the console during the search. Defaults to \code{FALSE}. }
}

} % END DETAILS

\value{A list with the following components by default. See the options of \code{control} for how to obtain a graphical representation of the derivation or how to benchmark the search. 

\describe{
  \item{\code{identifiable}}{A logical value that attains the value \code{TRUE} is the target quantity is identifiable and \code{FALSE} otherwise. }
  \item{\code{formula}}{A character string describing an identifiable query or an empty character vector for an unidentifiable effect. }
}

} % END VALUE

\references{

Bareinboim E., Pearl J. 2014 Transportability from Multiple Environments with Limited Experiments: Completeness Results. In \emph{Proceedings of the 27th Annual Conference on Neural Information Processing Systems}, 280--288.

Bareinboim E., Tian J. 2015 Recovering Causal Effects From Selection Bias. In \emph{Proceedings of the 29th AAAI Conference on Artificial Intelligence}, 3475--3481.

Mohan K., Pearl J., Tian, J. 2013 Graphical Models for Inference with Missing Data. In \emph{Advances in Neural Information Systems}, \bold{26}, 1277--1285.

Shpitser I., Pearl J. 2006 Identification of Joint Interventional Distributions in Recursive semi-Markovian Causal Models. In \emph{Proceedings of the 21st National Conference on Artificial Intelligence}, \bold{2}, 1219--1226. 

Shpitser I., Pearl J. 2006 Identification of Conditional Interventional Distributions. In \emph{Proceedings of the 22nd Conference on Uncertainty in Artificial Intelligence}, 427--444. 

} % END REFERENCES

\author{Santtu Tikka}

\examples{
library(dosearch)

# Multiple input distributions (both observational and interventional)

data1 <- "
  p(z_2,x_2|do(x_1))
  p(z_1|x_2,do(x_1,y))
  p(x_1|w_1,do(x_2))
  p(y|z_1,z_2,x_1,do(x_2))
  p(w|y,x_1,do(x_2))
"

query1 <- "p(y,x_1|w,do(x_2))"

graph1 <- "
  x_1 -> z_2
  x_1 -> z_1
  x_2 -> z_1
  x_2 -> z_2
  z_1 -> y
  z_2 -> y
  x_1 -> w
  x_2 -> w
  z_1 -> w
  z_2 -> w
"

get_derivation(data1, query1, graph1)

# Selection bias

data2 <- "
  p(x,y,z_1,z_2|s)
  p(z_1,z_2)
"

query2 <- "p(y|do(x))"

graph2 <- "
  x   -> z_1
  z_1 -> z_2
  x   -> y
  y   -- z_2
  z_2 -> s
"

get_derivation(data2, query2, graph2, selection_bias = "s")

# Transportability

data3 <- "
  p(x,y,z_1,z_2)
  p(x,y,z_1|s_1,s_2,do(z_2))
  p(x,y,z_2|s_3,do(z_1))
"

query3 <- "p(y|do(x))"

graph3 <- "
  z_1 -> x
  x   -> z_2
  z_2 -> y
  z_1 -- x
  z_1 -- z_2
  z_1 -- y
  s_1 -> z_1
  s_2 -> z_2
  s_3 -> y
"

get_derivation(data3, query3, graph3, transportability = "s_1, s_2, s_3")

# Missing data

data4 <- "
  p(x*,y*,z*,m_x,m_y,m_z)
"

query4 <- "p(x,y,z)"

graph4 <- "
  z -> x
  x -> y
  x -> m_z
  y -> m_z
  y -> m_x
  z -- y
"

get_derivation(data4, query4, graph4, missing_data = "m_x : x, m_y : y, m_z : z")

} % END EXAMPLES