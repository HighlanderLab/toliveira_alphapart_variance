## Functions

RecSysMaleFemale <- function(datafile, popname, gen, type = "Non-Selected") {
  datafile = rbind(datafile,
                   tibble(generation = gen,
                          ind        = popname@id,
                          father     = popname@father,
                          mother     = popname@mother,
                          sex        = popname@sex,
                          pheno      = popname@pheno[,1],
                          tbv        = popname@gv[,1],
                          type       = type))
}


meanP.2 = function(pop, na.rm = FALSE){   
  colMeans(pop@pheno, na.rm = na.rm)
}

PullSumm = function(datafile, popname, type = "Non-Selected"){
  gePa = genParam(popname)
  datafile = rbind(datafile,
                   tibble(generation = generation,
                          sex        = unique(popname@sex),
                          MeanGeno   = gePa$mu,
                          MeanA      = mean(gePa$gv_a),
                          VarA       = gePa$varA,
                          GenicVA    = gePa$genicVarA,
                          type       = type))
  return(datafile)
}

