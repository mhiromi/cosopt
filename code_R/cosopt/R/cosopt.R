`cosopt` <- function(data,sigma,timepoints,plotting=TRUE)
{
  library.dynam("cosopt_rlib")
  if(missing(data) || missing(sigma) || missing(timepoints)){
    stop(message="Please put a correct data.")
  }
  if(length(data) < 2){
    stop(message="Input time course data hase more than 2 points.")
  }
  if(length(data) != length(timepoints)){
    stop(message="Please equal number of timecourse column and number of timepoints column.")
  }
  if(length(data) != length(sigma)){
    stop(message="Please equal dimension of timecourse data and sigma data.")
  }
  if(!is.vector(data) || !is.vector(sigma) || !is.vector(timepoints)){
    stop(message="Each data must be given by vector.")
  }

  ans <- .C("cosopt",as.double(data),as.double(sigma),as.double(timepoints),as.integer(1),as.integer(length(data)),as.integer(plotting),ans=double(2),PACKAGE="cosopt")$ans
  if(plotting == TRUE && ans[1] != 0){
    freq <- ans[1]
    phase <- ans[2]
    test_cosine <- function(x) { cos(2*pi*freq*x+(2*pi*phase/100)) }
    plot(test_cosine,timepoints[1],timepoints[length(timepoints)])
    par(new=TRUE)
    plot(timepoints,data)
  }
}
