is.string = function(possible_string)
{
  return(is.character(possible_string) & length(possible_string) == 1)
}
  
#' ClonalOrigin.
#'
#' @param tree_filename The file containing the clonal frame tree.
#' @param data_filename The file containing the aligned sequence data.
#' @param output_filename The file to which the output will be written.
#' @param w Sets the number of pre burn-in iterations (default is 100000).
#' @param x Sets the number of burn-in iterations (default is 100000).
#' @param y Sets the number of iterations after burn-in (default is 100000).
#' @param z Sets the number of iterations between samples (default is 100).
#' @param T Sets the value of theta. Write "sNUM" in a string to set a per-site rate of NUM.
#' @param R Sets the value of rho. Write "sNUM" in a string to set a per-site rate of NUM.
#' @param D Sets the value of delta.
#' @param s Use given seed to initiate random number generator.
#' @param S Takes a string. If in the format "NUM,SEED", runs on a subset of NUM regions determined by seed SEED; if in the format "NUM/NUM/../NUM", runs on a specified region(s) given by each NUM.
#' @param r Perform r tempered steps between topological updates (default:0).
#' @param t Tempered at "temperature" t for topological updates (default:1).
#' @param U If TRUE, Start from UPGMA tree, rather than the default random tree.
#' @param G Greedily compute the "best fit" tree, given the recombination observed on the current tree.  If the input is negative and a previous run is provided, the tree is calculated from all observed values. If the input is positive, a "greedy move" is performed with this weight (see argument a). Note that this is NOT an MCMC move and causes bias.
#' @param a A vector or comma separated string setting the ELEVEN (real valued) move weightings to the given vector. The weightings need not sum to 1, but must be in the following order: MoveRho (ignored if not needed); MoveDelta (ignored if not needed); MoveTheta (ignored if not needed); MoveRemEdge; MoveAddEdge; MoveSiteChange; MoveTimeChange; MoveEdgeChange; MoveAgeClonal; MoveScaleTree; MoveRegraftClonal. 
#' @param m A vector or comma separated string setting THREE numbers required for MAIS implementation: gamma Power used in the annealing computation, gamma>1 means steps for jumping up concentrate near start, gamma=1 means steps are equidistant; T Number of steps for AIS, T=1 corresponds to no intermediate steps; N Number of replications for MAIS, instead of using N=1 (which corresponds to pure AIS-RJMCMC) use N=0.
#' @param i A vector or comma separated string setting the SIX parameters for creating random Recombination Trees under the inference model. The parameters are: N	(integer)	The number of sequences in the sample (default 10); n_B	(integer)	The number of block boundaries in the sample (default 8); l_B	(integer)	The length of each block: L=n_B * l_B (default 500); delta (real) The average length of imports (default 500.0); theta (real) The mutation rate NOT per site (default 100.0); rho	(real) The recombination rate NOT per site (default 50.0).
#' @param f If TRUE, forbid topology changes, (allowing updates of coalescence times).
#' @param v If TRUE, use verbose mode.
#' @param h If TRUE, display a help message.
#' @param V If TRUE, print ClonalOrigin Version info.
#' @return Writes output to a file specified by output_filename.
#' @export
clonal_origin = function(tree_filename,
                          data_filename,
                          output_filename,
                          w=NA,
                          x=NA,
                          y=NA,
                          z=NA,
                          T=NA,
                          R=NA,
                          D=NA,
                          s=NA,
                          S=NA,
                          r=NA,
                          t=NA,
                          U=NA,
                          G=NA,
                          a=NA,
                          m=NA,
                          i=NA,
                          f=NA,
                          v=NA,
                          h=NA,
                          V=NA)
{
  if (!file.exists(tree_filename))
    stop("Tree file does not exist.")
  
  if (!file.exists(data_filename)) 
    stop("Data file does not exist.")
  
  argv_list = list()
  
  if (!is.na(w))
  {
    if (is.numeric(w) && (w>=0) )
      argv_list[["w"]] = paste("-w",round(w),sep=" ")
    else
      stop("Number of pre-burn-in iterations must be non-negative.")
  }
  
  if (!is.na(x))
  {
    if (is.numeric(x) && (x>=0) )
      argv_list[["x"]] = paste("-x",round(x),sep=" ")
    else
      stop("Number of burn-in iterations must be non-negative.")
  }
  
  if (!is.na(y))
  {
    if (is.numeric(y) && (y>=0) )
      argv_list[["y"]] = paste("-y",round(y),sep=" ")
    else
      stop("Number of iterations after burn-in must be non-negative.")
  }
  
  if (!is.na(z))
  {
    if (is.numeric(z) && (z>=0) )
      argv_list[["z"]] = paste("-z",round(z),sep=" ")
    else
      stop("Number of iterations between samples must be non-negative.")
  }
  
  if (!is.na(T))
  {
    if (is.numeric(T)) {
      if (T>=0)
      {
        argv_list[["T"]] = paste("-T",T,sep=" ")
      }
      else
      {
        stop("theta must be non-negative.")
      }
    } else if (is.string(T)) {
      if (substr(T,1,1)=="s")
      {
        nch = nchar(T)
        if (nch>1)
        {
          if (as.double(substr(T,2,nch))>=0)
          {
            argv_list[["T"]] = paste("-T",T,sep=" ")
          }
          else
          {
            stop("theta must be non-negative.")
          }
        }
        else
        {
          stop("No theta given.")
        }
      }
      else
      {
        if (as.double(T)>=0)
        {
          argv_list[["T"]] = paste("-T",T,sep=" ")
        }
        else
        {
          stop("theta must be non-negative.")
        }
      }
    }
    else {
      stop("Argument must be a float or a string.")
    }
  }
  
  if (!is.na(R))
  {
    if (is.numeric(R)) {
      if (R>=0)
      {
        argv_list[["R"]] = paste("-R",R,sep=" ")
      }
      else
      {
        stop("rho must be non-negative.")
      }
    } else if (is.string(R)) {
      if (substr(R,1,1)=="s")
      {
        nch = nchar(R)
        if (nch>1)
        {
          if (as.double(substr(R,2,nch))>=0)
          {
            argv_list[["R"]] = paste("-R",R,sep=" ")
          }
          else
          {
            stop("rho must be non-negative.")
          }
        }
        else
        {
          stop("No rho given.")
        }
      }
      else
      {
        if (as.double(R)>=0)
        {
          argv_list[["R"]] = paste("-R",R,sep=" ")
        }
        else
        {
          stop("rho must be non-negative.")
        }
      }
    }
    else {
      stop("Argument must be a float or a string.")
    }
  }
  
  if (!is.na(D))
  {
    if (is.numeric(D) && (D>=0) )
      argv_list[["D"]] = paste("-D",D,sep=" ")
    else
      stop("delta must be non-negative.")
  }
  
  if (!is.na(s))
  {
    if (is.numeric(s) && (s>=0) )
      argv_list[["s"]] = paste("-s",round(s),sep=" ")
    else
      stop("Seed must be non-negative.")
  }
  
  if (!is.na(S))
  {
    if (is.string(S))
    {
      if ( grepl(",",S) || grepl("/",S) )
      {
        argv_list[["S"]] = paste("-S",S,sep=" ")
      }
      else
      {
        stop("Invalid subset of regions.")
      }
    }
    else
    {
      stop("Subset of regions must be specified using a string.")
    }
  }
  
  if (!is.na(r))
  {
    if (is.numeric(r) && (r>=0) )
      argv_list[["r"]] = paste("-r",round(r),sep=" ")
    else
      stop("Number of tempered steps between topological updates must be non-negative.")
  }
  
  if (!is.na(t))
  {
    if (is.numeric(t) && (t>=0) )
      argv_list[["t"]] = paste("-t",r,sep=" ")
    else
      stop("Temperature for topological updates must be non-negative.")
  }
  
  if (!is.na(U))
  {
    if (U==TRUE)
      argv_list[["U"]] = "-U"
  }
  
  if (!is.na(G))
  {
    if (is.numeric(G)>0)
      argv_list[["G"]] = paste("-G",G,sep=" ")
    else
      stop("Temperature for topological updates must be positive.")
  }
  
  if (!is.na(a))
  {
    if (is.string(a))
    {
      if (grepl(",",a))
      {
        split_a = strsplit(a,",")
        if (length(split_a[[1]])!=11)
          stop("Require 11 move weightings.")
        if ( (as.numeric(split_a[[1]][1])>=0) && (as.numeric(split_a[[1]][2])>=0) && (as.numeric(split_a[[1]][3])>=0) && (as.numeric(split_a[[1]][4])>=0) && (as.numeric(split_a[[1]][5])>=0) && (as.numeric(split_a[[1]][6])>=0) && (as.numeric(split_a[[1]][7])>=0) && (as.numeric(split_a[[1]][8])>=0) && (as.numeric(split_a[[1]][9])>=0) && (as.numeric(split_a[[1]][10])>=0) && (as.numeric(split_a[[1]][11])>=0) )
          argv_list[["a"]] = paste("-a",a,sep=" ")
        else
          stop("Move weightings must be non-negative.")
      }
      else
      {
        stop("Invalid move weightings.")
      }
    }
    else if (is.vector(a))
    {
      if (length(a)!=11)
        stop("Require 11 move weightings.")
      if ( (a[1]>=0) && (a[2]>=0) && (a[3]>=0) && (a[4]>=0) && (a[5]>=0) && (a[6]>=0) && (a[7]>=0) && (a[8]>=0) && (a[9]>=0) && (a[10]>=0) && (a[11]>=0) )
        argv_list[["a"]] = paste("-a",paste(a,collapse = ","),sep=" ")
      else
        stop("Move weightings must be non-negative.")
    }
    else
    {
      stop("Move weightings must be specified using a string or a vector.")
    }
  }
  
  if (!is.na(m))
  {
    if (is.string(m))
    {
      if (grepl(",",m))
      {
        split_m = strsplit(m,",")
        if (length(split_m[[1]])!=3)
          stop("Require 3 MAIS parameters.")
        if ( (as.numeric(split_m[[1]][1])>0) && (as.numeric(split_m[[1]][2])>=1) && (as.numeric(split_m[[1]][3])>=0) )
          argv_list[["m"]] = paste("-m",m,sep=" ")
        else
          stop("Invalid MAIS parameters.")
      }
      else
      {
        stop("Invalid move weightings.")
      }
    }
    else if (is.vector(m))
    {
      if (length(m)!=3)
        stop("Require 3 move weightings.")
      if ( (m[1]>0) && (m[2]>=1) && (m[3]>=0) )
        argv_list[["m"]] = paste("-m",paste(m,collapse = ","),sep=" ")
      else
        stop("Invalid MAIS parameters.")
    }
    else
    {
      stop("MAIS parameters must be specified using a string or a vector.")
    }
  }
  
  if (!is.na(i))
  {
    if (is.string(i))
    {
      if (grepl(",",i))
      {
        split_i = strsplit(i,",")
        if (length(split_i[[1]])!=6)
          stop("Require 6 parameters for creating recombination trees.")
        if ( (as.numeric(split_i[[1]][1])>0) && (as.numeric(split_i[[1]][2])>=0) && (as.numeric(split_i[[1]][3])>=1) && (as.numeric(split_i[[1]][4])>=0) && (as.numeric(split_i[[1]][5])>=0) && (as.numeric(split_i[[1]][6])>=0) )
          argv_list[["i"]] = paste("-i",i,sep=" ")
        else
          stop("Invalid parameters for creating recombination trees.")
      }
      else
      {
        stop("Invalid parameters for creating recombination trees.")
      }
    }
    else if (is.vector(i))
    {
      if (length(i)!=6)
        stop("Require 6 parameters for creating recombination trees.")
      if ( (i[1]>0) && (i[2]>=0) && (i[3]>=1) && (i[4]>=0) && (i[5]>=0) && (i[6]>=0) )
        argv_list[["i"]] = paste("-i",paste(i,collapse = ","),sep=" ")
      else
        stop("Invalid parameters for creating recombination trees.")
    }
    else
    {
      stop("Parameters for creating recombination trees must be specified using a string or a vector.")
    }
  }
  
  if (!is.na(f))
  {
    if (f==TRUE)
      argv_list[["f"]] = "-f"
  }
  
  if (!is.na(v))
  {
    if (v==TRUE)
      argv_list[["v"]] = "-v"
  }
  
  if (!is.na(h))
  {
    if (h==TRUE)
      argv_list[["h"]] = "-h"
  }
  
  if (!is.na(V))
  {
    if (V==TRUE)
      argv_list[["V"]] = "-V"
  }
  
  # if (!is.na(L))
  # {
  #   if (L==TRUE)
  #     argv_list[[L]] = "-L"
  # }
  # 
  # if (!is.na(C))
  # {
  #   if (C==TRUE)
  #     argv_list[[C]] = "-C"
  # }
  
  argv_list[[length(argv_list)+1]] = tree_filename
  argv_list[[length(argv_list)+1]] = data_filename
  argv_list[[length(argv_list)+1]] = output_filename
  
  run_clonal_origin(argv_list)
}