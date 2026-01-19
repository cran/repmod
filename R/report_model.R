#' Report from linear model
#'
#' @description Creates a report table from a linear model.
#' @param x A linear model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats confint getCall
#' @export
report.lm<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  sx <- summary(x)
  ci <- confint(x)
  obj <- list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx$coefficients[,2], lwr.int=ci[,1],
              upper.int=ci[,2], pvalues=sx$coefficients[,4], r.squared=sx$r.squared, adj.r.squared=sx$adj.r.squared, sigma=sx$sigma)
  output<-rbind(cbind(round(obj$coefficients,digits),round(obj$se,digits),
                      round(obj$lwr.int,digits),round(obj$upper.int, digits), round(obj$pvalues,digitspvals)),
                c(round(obj$r.squared,digits+1),rep("",4)),
                c(round(obj$adj.r.squared,digits+1),rep("",4)),
                c(round(obj$sigma, digits+1), rep("", 4)))
  colnames(output)<-c('Estimate','Std. Error','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$coefficients)[1]+c(1,2,3)]<-c('R Squared','Adj.R Squared', "Sigma")
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Report from generalized linear model
#'
#' @description Creates a report table from a generalized linear model.
#' @param x A generalized linear model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.glm<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  compute.exp<-x$family$link %in% c("logit", "log")
  sx<-summary(x)
  ci<-confint(x)
  obj <- list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx$coefficients[,2], lwr.int=ci[,1],
              upper.int=ci[,2], pvalues=sx$coefficients[,4], aic=sx$aic)
  if(compute.exp){
    obj$exp.coef <- exp(sx$coefficients[,1])
    obj$exp.lwr.int <- exp(ci[,1])
    obj$exp.upper.int <- exp(ci[,2])
  }
  output<-rbind(cbind(round(obj$coefficients,digits), round(obj$se, digits),
                      if(compute.exp) {
                        cbind(round(obj$exp.coef,digits), round(obj$exp.lwr.int, digits),
                              round(obj$exp.upper.int, digits))
                      } else{
                        cbind(round(obj$lwr.int, digits), round(obj$upper.int, digits))
                      }
                      , round(obj$pvalues,digitspvals)),
                c(round(obj$aic,digits),rep("",ifelse(compute.exp, 5, 4))))
  colnames(output)<-c('Estimate','Std. Error',if(compute.exp) {'exp(Estimate)'},'Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$coefficients)[1]+1]<-c('AIC')
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Report from cox regression model
#'
#' @description Creates a report table from a cox model.
#' @param x A cox model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats AIC getCall
#' @export
report.coxph<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  sx<-summary(x)
  obj <- list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx$coefficients[,3], hr=sx$conf.int[,1],
              lwr.int=sx$conf.int[,3], upper.int=sx$conf.int[,4], pvalues=sx$coefficients[,5], aic=AIC(x))
  output<-rbind(cbind(round(obj$coefficients,digits), round(obj$se, digits), round(obj$hr, digits),
                      round(obj$lwr.int,digits), round(obj$upper.int, digits),
                      round(obj$pvalues,digitspvals)),c(round(obj$aic,digits),rep('',5)))
  colnames(output)<-c('Estimate','Std. Error','HR','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(output)[1]]<-c('AIC')
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Report from linear mixed model with pvalues
#'
#' @description Creates a report table from a linear mixed model.
#' @param x A linear mixed model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @importFrom methods loadMethod
#' @export
report.merModLmerTest<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  #loadNamespace("lmerTest")
  loadMethod("summary", "summary.lmerModLmerTest", envir="lmerTest")
  sx <- summary(x)
  #unloadNamespace("lmerTest")
  cor <- as.data.frame(lme4::VarCorr(x))
  ci <- confint(x)
  #cor[dim(cor)[1],2]<-'Residual'
  obj<- list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx$coefficients[,2], lwr.int=ci[,1][-c(1:dim(as.data.frame(lme4::VarCorr(x)))[1])],
             upper.int=ci[,2][-c(1:dim(as.data.frame(lme4::VarCorr(x)))[1])],
             pvalues=tryCatch(sx$coefficients[,5], error=function(x) NA), aic=AIC(x),
             random=cor[c(is.na(cor$var2)),c(5)])
  output<-rbind(rbind(cbind(round(obj$coefficients,digits),round(obj$se,digits),
                            round(obj$lwr.int,digits), round(obj$upper.int, digits),
                            round(obj$pvalues, digits)),
                      c(round(obj$aic,digits-1),rep("",4))),
                matrix(c(round(obj$random,digits),rep("",4*dim(cor[c(is.na(cor$var2)),])[1])),ncol=5,byrow=F))
  colnames(output)<-c('Estimate','Std. Error','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$coefficients)[1]+1]<-c('AIC')
  rownames(output)[rownames(output)==""]<-paste(c(rep('Sd ',length(cor[is.na(cor$var2), c(2)])-1), ""),
                                                na.omit(cor[is.na(cor$var2), c(1)]),
                                                c(na.omit(cor[is.na(cor$var2), c(2)]), ""),sep='')
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Report from linear mixed model
#'
#' @description Creates a report table from a linear mixed model.
#' @param x A linear mixed model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.lmerMod<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  x<-lmerTest::lmer(x@call,data=x@frame)
  report.merModLmerTest(x, file, type, digits, digitspvals, info=info, ...)
}

#' Report from generalized linear mixed model
#'
#' @description Creates a report table from a generalized linear mixed model.
#' @param x A generalized linear mixed model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.glmerMod<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  compute.exp<-x@resp$family$link %in% c("logit", "log")
  sx<-summary(x)
  cor<-as.data.frame(lme4::VarCorr(x))
  ci <- confint(x)
  obj<- list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx$coefficients[,2], lwr.int=ci[,1][-c(1:dim(as.data.frame(lme4::VarCorr(x)))[1])],
             upper.int=ci[,2][-c(1:dim(as.data.frame(lme4::VarCorr(x)))[1])],
             pvalues=sx$coefficients[,4], aic=AIC(x),
             random=cor[c(is.na(cor$var2)),c(5)])
  if(compute.exp){
    obj$exp.coef <- exp(obj$coefficients)
    obj$exp.lwr.int <- exp(obj$lwr.int)
    obj$exp.upper.int <- exp(obj$upper.int)
  }
  output<-rbind(rbind(cbind(round(obj$coefficients,digits),
                            round(obj$se, digits),
                            if(compute.exp) {
                              cbind(round(obj$exp.coef,digits), round(obj$exp.lwr.int, digits),
                                    round(obj$exp.upper.int, digits))
                            } else{
                              cbind(round(obj$lwr.int, digits), round(obj$upper.int, digits))
                            }
                            , round(obj$pvalues,digitspvals)),
                      c(round(obj$aic, digits), rep("", ifelse(compute.exp, 5, 4)))),
                matrix(c(round(obj$random,digits),rep("",ifelse(compute.exp, 5, 4)*dim(cor[c(is.na(cor$var2)),])[1])),ncol=ifelse(compute.exp, 6, 5),byrow=F))
  colnames(output)<-c('Estimate','Std. Error', if(compute.exp) 'exp(Estimate)','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$coefficients)[1]+1]<-c('AIC')
  rownames(output)[rownames(output)==""]<-paste(rep('Sd ',length(cor[is.na(cor$var2), c(2)])),
                                                na.omit(cor[is.na(cor$var2), c(1)]),
                                                na.omit(cor[is.na(cor$var2), c(2)]),sep='')
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}


#' Report from quantile mixed model
#'
#' @description Creates a report table from a quantile mixed model.
#' @param x A quantile model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.lqmm<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  sx<-summary(x, ...)
  obj<-list(coefficients=setNames(sx$tTable[,1], rownames(sx$tTable)), se=sx$tTable[,2], lwr.int=sx$tTable[,3], upper.int=sx$tTable[,4],
            pvalues=sx$tTable[,5], aic=sx$aic, random=round(VarCorr(x),2))
  output<-rbind(rbind(cbind(round(obj$coefficients,digits), round(obj$se,digits),
                            round(obj$lwr.int, digits), round(obj$upper.int, digits), round(obj$pvalues, digitspvals)),
                      c(round(obj$aic, digits),rep("",4))),
                matrix(c(round(obj$random,digits),rep("",4*length(obj$random))),ncol=5,byrow=F))
  colnames(output)<-c('Estimate','Std. Error','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$tTable)[1]+1]<-c('AIC')
  rownames(output)[(dim(sx$tTable)[1]+2):(((dim(sx$tTable)[1]+2)+length(obj$random))-1)]<-paste('Ran.Eff',names(VarCorr(x)),sep=' ')
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}


#' Report from ordinal model
#'
#' @description Creates a report table from an ordinal model.
#' @param x An ordinal model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.clm<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  compute.exp<-x$link %in% c("logit", "log")
  sx<-summary(x)
  ci<-confint(x)
  obj<-list(coefficients=sx$coefficients[,1][-c(1:length(x$alpha))], se=sx$coefficients[,2][-c(1:length(x$alpha))],
            lwr.int=ci[,1], upper.int=ci[,2], pvalues=sx$coefficients[,4][-c(1:length(x$alpha))], aic=AIC(x))
  if(compute.exp){
    obj$exp.coef <- exp(obj$coefficients)
    obj$exp.lwr.int <- exp(obj$lwr.int)
    obj$exp.upper.int <- exp(obj$upper.int)
  }
  obj$thresholds <- sx$coefficients[1:length(x$alpha),]
  output<-rbind(cbind(round(obj$coefficients,digits), round(obj$se,digits),
                      if(compute.exp) {
                        cbind(round(obj$exp.coef, digits), round(obj$exp.lwr.int,digits),
                              round(obj$exp.upper.int, digits))
                      } else{
                        cbind(round(obj$lwr.int, digits), round(obj$upper.int, digits))
                      }
                      , round(obj$pvalues,digitspvals)), c(round(obj$aic,digits),rep("", ifelse(compute.exp, 5, 4))))
  colnames(output)<-c('Estimate','Std. Error',if(compute.exp) 'exp(Estimate)','Lower 95%','Upper 95%','P-value')
  rownames(output)[length(rownames(output))]<-c('AIC')
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Report from ordinal mixed model
#'
#' @description Creates a report table from an ordinal mixed model.
#' @param x An ordinal model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.clmm<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  compute.exp<-x$link %in% c("logit", "log")
  sx<-summary(x)
  ci<-confint(x)
  obj <- list(coefficients=sx$coefficients[,1][-c(1:length(x$alpha))], se=sx$coefficients[,2][-c(1:length(x$alpha))],
              lwr.int=ci[,1][-c(1:length(x$alpha))], upper.int=ci[,2][-c(1:length(x$alpha))], pvalues=sx$coefficients[,4][-c(1:length(x$alpha))], aic=AIC(x),
              random=VarCorr(x))
  if(compute.exp){
    obj$exp.coef <- exp(obj$coefficients)
    obj$exp.lwr.int <- exp(obj$lwr.int)
    obj$exp.upper.int <- exp(obj$upper.int)
  }
  obj$thresholds <- sx$coefficients[1:length(x$alpha),]
  output<-rbind(rbind(cbind(round(obj$coefficients,digits), round(obj$se, digits),
                            if(compute.exp) {
                              cbind(round(obj$exp.coef,digits), round(obj$exp.lwr.int, digits),
                                    round(obj$exp.upper.int, digits))
                            } else{
                              cbind(round(obj$lwr.int, digits), round(obj$upper.int, digits))
                            }
                            , round(obj$pvalues,digitspvals)), c(round(obj$aic,digits),rep("",ifelse(compute.exp, 5, 4)))),
                matrix(c(round(rapply(obj$random, function(x) sqrt(diag(x))),digits),
                         rep("",ifelse(compute.exp, 5, 4)*length(rapply(obj$random, function(x) sqrt(diag(x)))))),ncol=ifelse(compute.exp, 6, 5),byrow=F))
  colnames(output)<-c('Estimate','Std. Error', if(compute.exp) 'exp(Estimate)','Lower 95%','Upper 95%','P-value')
  rownames(output)[length(x$beta)+1]<-c('AIC')
  rownames(output)[rownames(output)==""]<-names(rapply(VarCorr(x), function(x) sqrt(diag(x))))
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}


#' Report from quantile regression model
#'
#' @description Creates a report table from a quantile regression model.
#' @param x A quantreg model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.rq<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  sx<-summary(x, se="rank")
  sx2<-summary(x, covariance=TRUE)
  obj<-list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx2$coefficients[,2], lwr.int=sx$coefficients[,2],
            upper.int=sx$coefficients[,3], pvalues=sx2$coefficients[,4], aic=AIC(x))
  output<-rbind(cbind(round(obj$coefficients,digits),round(obj$se,digits),
                      round(obj$lwr.int,digits), round(obj$upper.int, digits),
                      round(obj$pvalues,digitspvals)),
                c(round(obj$aic, digits),rep("",4)))
  colnames(output)<-c('Estimate','Std. Error','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$coefficients)[1]+1]<-'AIC'
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}


#' Report from beta regression model
#'
#' @description Creates a report table from a beta regression model.
#' @param x A betareg model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.betareg<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  compute.exp<-x$link$mean$name %in% c("logit", "log")
  sx<-summary(x)
  ci<-confint(x)
  obj<-list(coefficients=setNames(sx$coefficients$mean[,1], rownames(sx$coefficients$mean)), se=sx$coefficients$mean[,2], lwr.int=ci[,1][-dim(ci)[1]],
            upper.int=ci[,2][-dim(ci)[1]], pvalues=sx$coefficients$mean[,4], aic=AIC(x), pseudo.r=sx$pseudo.r.squared,
            phi=sx$coefficients$precision)
  if(compute.exp){
    obj$exp.coef <- exp(obj$coefficients)
    obj$exp.lwr.int <- exp(obj$lwr.int)
    obj$exp.upper.int <- exp(obj$upper.int)
  }
  output<-rbind(cbind(round(obj$coefficients,digits),round(obj$se,digits),
                      if(compute.exp) {
                        cbind(round(obj$exp.coef,digits), round(obj$exp.lwr.int, digits),
                              round(obj$exp.upper.int, digits))
                      } else{
                        cbind(round(obj$lwr.int, digits), round(obj$upper.int, digits))
                      }
                      , round(obj$pvalues, digitspvals)), c(round(obj$phi[1], digits+1),rep("", ifelse(compute.exp, 5, 4))),
                c(round(obj$pseudo.r, digits), rep("", ifelse(compute.exp, 5, 4))))
  colnames(output)<-c('Estimate','Std. Error',if(compute.exp) 'exp(Estimate)', 'Lower 95%','Upper 95%','P-value')
  rownames(output)[c(dim(sx$coefficients$mean)[1]+1, dim(sx$coefficients$mean)[1]+2)]<-c("phi", "Pseudo R-squared")
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Report models from brms package
#'
#' @description Creates a report table from model fitted by brms.
#' @param x A brms model object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param prob Probability level for credible intervals
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @importFrom bayestestR pd
#' @export
report.brmsfit <- function(x, file=NULL, type="word", digits=3, info=TRUE, print=TRUE, exclude=NULL, prob=0.95, ...){
  compute.exp<-x$family$link %in% c("logit", "log")
  sx<-summary(x, prob=prob)
  WC<-eval(parse(text="brms::WAIC(x)"))
  random<-tryCatch(do.call(rbind, sx$random), error=function(e) NA)
  pd <- bayestestR::pd(x)
  #if(!any(is.na(random))) rownames(random)<-paste(rownames(random),rep(names(sx$random), sapply(sx$random, nrow)), sep=" ")
  obj<-list(coefficients=setNames(sx$fixed[,1], rownames(sx$fixed)), se=sx$fixed[,2], lwr.int=sx$fixed[,3],
            upper.int=sx$fixed[,4], pd=pd$pd, random=random, WAIC=setNames(c(WC$estimates[3,1], WC$estimates[3,2]), c("WAIC", "WAIC SE")), Eff.Sample_min=min(c(sx$fixed[,6], sx$fixed[,7])), Rhat_max=round(max(sx$fixed[,5]),2))
  if(compute.exp){
    obj$exp.coef <- exp(obj$coefficients)
    obj$exp.lwr.int <- exp(obj$lwr.int)
    obj$exp.upper.int <- exp(obj$upper.int)
  }
  output<-rbind(cbind(round(obj$coefficients,digits),round(obj$se,digits),
                      if(compute.exp){
                        cbind(round(obj$exp.coef, digits), round(obj$exp.lwr.int, digits),
                              round(obj$exp.upper.int, digits))
                      } else{
                        cbind(round(obj$lwr.int,digits), round(obj$upper.int, digits))
                      }, round(pd$pd, digits)), if(!any(is.na(random))) {
                        if(compute.exp){
                          as.matrix(cbind(round(random[,1:2, drop=FALSE], digits), "-", round(random[,3:4, drop=FALSE], digits), ""))
                        } else {
                          as.matrix(cbind(round(random[,1:2, drop=FALSE], digits), round(random[,3:4, drop=FALSE], digits), ""))
                        }
                      },
                c(round(WC$estimates[3,1], digits), round(WC$estimates[3,2], digits), rep("", ifelse(compute.exp, 3, 2)), ""))
  rownames(output)[dim(output)[1]]<-"WAIC"
  colnames(output)<-c('Estimate','Std. Error',if(compute.exp) 'exp(Estimate)', paste('Lower ', round(prob*100), '%', sep=""),paste('Upper ', round(prob*100), '%', sep=""), "pd")
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    suppressWarnings(make_table(output, file, type, info=info, ...))
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  if(obj$Rhat_max > 1.1) warning("Please diagnose your model, Rhat values greater than 1.1")
  class(obj) <- "reportmodel"
  invisible(obj)
}


#' Report models from glmnet package
#'
#' @description Creates a report table from models fitted by glmnet.
#' @param x A glmnet model object
#' @param s Value of lambda for estimating the coefficients
#' @param gamma Value of gamma for estimating the coefficients (only used in relaxed fits)
#' @param drop.zero Should zero coefficients be dropped?
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats coef getCall
#' @export
report.glmnet<-function(x, s, gamma=1, drop.zero=TRUE, file=NULL, type="word", digits=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  compute.exp<- any(grepl("binomial|cox", x$call))
  if("relaxed" %in% class(x)){
    coefs <- coef(x, s=s, gamma=gamma)
  } else {
    coefs <- coef(x, s=s)
  }
  if("multnet" %in% class(x)){
    temp_coef <- setNames(data.frame(as.matrix(do.call(cbind, coefs))), names(coefs))
    obj <- list(coefficients=temp_coef[apply(temp_coef, 1, function(x) any(x != 0)),], lwr.int=NA, upper.int=NA)
  } else {
    obj <- list(coefficients=as.numeric(coefs)[if(drop.zero) {as.numeric(coefs)!=0}], lwr.int=NA, upper.int=NA)
    names(obj$coefficients)<-rownames(coefs)[if(drop.zero) {as.numeric(coefs)!=0}]
  }
  if(compute.exp){
    obj$exp.coef <- exp(obj$coefficients)
  }
  obj$lambda <- s
  if("multnet" %in% class(x)){
    output <- rbind(round(obj$coefficients, digits), c(s, rep("", ifelse(compute.exp, max(ncol(obj$coefficients), 1), max(ncol(obj$coefficients)-1, 0)))))
    rownames(output)[nrow(output)] <- "lambda"
  } else {
    output<-rbind(cbind(round(obj$coefficients,digits),
                        if(compute.exp){
                          round(obj$exp.coef, digits)
                        }), cbind(s, rep("", ifelse(compute.exp, 1, 0))))
    colnames(output)<-c('Estimate', if(compute.exp) 'exp(Estimate)')
    rownames(output)<-c(names(obj$coefficients), "lambda")
  }
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}


#' Report from robust linear model (rlm)
#'
#' @description Creates a report table from a robust linear model.
#' @param x A rlm object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.rlm<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  sx <- summary(x, method = "XtWX")
  ci <- rob.ci(x, ...)
  obj <- list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx$coefficients[,2], lwr.int=ci[,1],
              upper.int=ci[,2], pvalues=rob.pvals(x), AIC=AIC(x))
  output<-rbind(cbind(round(obj$coefficients,digits),round(obj$se,digits),
                      round(obj$lwr.int,digits),round(obj$upper.int, digits), round(obj$pvalues,digitspvals)),
                c(round(obj$AIC,digits+1),rep("",4)))
  colnames(output)<-c('Estimate','Std. Error','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$coefficients)[1]+1]<-'AIC'
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}


#' Report from generalized linear mixed model from ADMB
#'
#' @description Creates a report table from a glmmadmb model.
#' @param x A generalized linear mixed model object (glmmabmb)
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats getCall
#' @export
report.glmmadmb<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  compute.exp<-x$link %in% c("logit", "log")
  sx<-summary(x)
  cor<-sqrt(cbind(unlist(lapply(x$S, function(x) diag(x)))))
  ci <- confint(x)
  obj<- list(coefficients=setNames(sx$coefficients[,1], rownames(sx$coefficients)), se=sx$coefficients[,2], lwr.int=ci[,1],
             upper.int=ci[,2],
             pvalues=sx$coefficients[,4], aic=AIC(x),
             random=cor)
  if(compute.exp){
    obj$exp.coef <- exp(obj$coefficients)
    obj$exp.lwr.int <- exp(obj$lwr.int)
    obj$exp.upper.int <- exp(obj$upper.int)
  }
  output<-rbind(rbind(cbind(round(obj$coefficients,digits),
                            round(obj$se, digits),
                            if(compute.exp) {
                              cbind(round(obj$exp.coef,digits), round(obj$exp.lwr.int, digits),
                                    round(obj$exp.upper.int, digits))
                            } else{
                              cbind(round(obj$lwr.int, digits), round(obj$upper.int, digits))
                            }
                            , round(obj$pvalues,digitspvals)),
                      c(round(obj$aic, digits), rep("", ifelse(compute.exp, 5, 4)))),
                matrix(c(round(obj$random,digits),rep("",ifelse(compute.exp, 5, 4)*dim(cor)[1])),ncol=ifelse(compute.exp, 6, 5),byrow=F))
  colnames(output)<-c('Estimate','Std. Error', if(compute.exp) 'exp(Estimate)','Lower 95%','Upper 95%','P-value')
  rownames(output)[dim(sx$coefficients)[1]+1]<-c('AIC')
  rownames(output)[rownames(output)==""]<-paste(rep('Sd ',dim(cor)[1]),
                                                rownames(cor),sep='')
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Report from generalized least squares model
#'
#' @description Creates a report table from a generalized least squares model.
#' @param x A gls object
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimals
#' @param digitspvals Number of decimals for p-values
#' @param info If TRUE, include call in the exported table
#' @param print Should the report table be printed on screen?
#' @param exclude Vector with rownames to remove from output
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @importFrom stats confint getCall
#' @export
report.gls<-function(x, file=NULL, type="word", digits=3, digitspvals=3, info=TRUE, print=TRUE, exclude=NULL, ...){
  sx <- summary(x)
  ci <- confint(x)
  obj <- list(coefficients=setNames(sx$tTable[,1], rownames(sx$tTable)), se=sx$tTable[,2], lwr.int=ci[,1],
              upper.int=ci[,2], pvalues=sx$tTable[,4], corStruct=coef(sx$modelStruct$corStruct, unconstrained = FALSE),
              varStruct=coef(sx$modelStruct$varStruct, unconstrained = FALSE), aic=sx$AIC, sigma=sx$sigma)
  output<-rbind(cbind(round(obj$coefficients,digits),round(obj$se,digits),
                      round(obj$lwr.int,digits),round(obj$upper.int, digits), round(obj$pvalues,digitspvals)),
                if(!is.null(sx$modelStruct$corStruct)) matrix(c(round(coef(sx$modelStruct$corStruct, unconstrained = FALSE), digits), rep("", 4*length(coef(sx$modelStruct$corStruct, unconstrained = FALSE)))), ncol=5),
                if(!is.null(sx$modelStruct$varStruc)) matrix(c(round(coef(sx$modelStruct$varStruct, unconstrained = FALSE), digits), rep("", 4*length(coef(sx$modelStruct$varStruct, unconstrained = FALSE)))), ncol=5),
                c(round(obj$aic, digits),rep("", 4)),
                c(round(obj$sigma, digits), rep("", 4)))
  colnames(output)<-c('Estimate','Std. Error','Lower 95%','Upper 95%','P-value')
  rownames(output)[-(1:dim(sx$tTable)[1])] <- c(if(!is.null(sx$modelStruct$corStruct)) names(coef(sx$modelStruct$corStruct, unconstrained = FALSE)),
                                                if(!is.null(sx$modelStruct$varStruct)) names(coef(sx$modelStruct$varStruct, unconstrained = FALSE)),
                                                'AIC', "Sigma")
  output[,"P-value"][output[,"P-value"]=="0"]<-"<0.001"
  output <- output[!rownames(output) %in% exclude,]
  if(!is.null(file)){
    info <- if(info) deparse1(getCall(x)) else NULL
    make_table(output, file, type, info=info, ...)
  }
  obj$output <- data.frame(output, check.names=FALSE, stringsAsFactors=FALSE)
  if(print) print(obj$output, row.names=TRUE, right=TRUE)
  class(obj) <- "reportmodel"
  invisible(obj)
}

#' Function to compute p-values for robust linear regression models
#'
#' @description Estimates p-values for rlm models.
#' @param x A rlm object
#' @return A vector of p-values
#' @importFrom stats pf
#' @export
rob.pvals <- function(x){
  coefs <- x$coef
  sx<-summary(x, method = "XtWX")
  covs<-diag(sx$cov.unscaled)*sx$stddev^2
  statistics<-sapply(1:length(coefs), function(x) sum(coefs[x]*solve(covs[x], coefs[x])))
  pf(statistics, 1, sx$df[2], lower.tail = FALSE)
}

#' Function to compute bootstrap confidence intervals for robust linear regression models
#'
#' @description Estimates confidence intervals for rlm models.
#' @param x A rlm object
#' @param level Confidence level for the interval
#' @param maxit Maximum number of iterations per fit
#' @param R Number of boostrap samples
#' @return A matrix with bootstrap confidence intervals for each variable in the model
#' @importFrom stats formula
#' @export
rob.ci <- function(x, level=0.95, maxit=200, R=2000){
  coefb <- function(object, data, indices){      #Function to extract R2
    d <- data[indices,]
    object$call$data <- as.name("d")
    object$call$maxit <- maxit
    fit <- eval(object$call)
    return(coef(fit))
  }
  results <- boot::boot(data=x$model, statistic=coefb, R=R, object=x)
  t(sapply(lapply(1:length(x$coefficients), function(x) boot::boot.ci(results, conf=level, type="bca", index=x)), function(x) x$bca[4:5]))
}

#' Export a table to word
#'
#' @description Exports a table to Word.
#' @param x A data.frame object
#' @param file Name of the file
#' @param info Footer for the table
#' @param use.rownames Should row names be added to the output?
#' @return Creates a word file with the table
#' @importFrom stats median na.omit quantile sd setNames
#' @export
make_word_table <- function(x, file, info=NULL, use.rownames=TRUE){
  mydoc <- officer::read_docx()
  dataf <- data.frame(x)
  if(use.rownames) dataf <- data.frame(Variables=rownames(dataf), dataf)
  MyFlexTable <- flextable::flextable(dataf)
  MyFlexTable <- flextable::bold(flextable::fontsize(flextable::font(MyFlexTable,fontname='Calibri',part='all'),
                                                     size = 12, part = "all"), bold = TRUE, part = "header")
  MyFlexTable <- set_noms(MyFlexTable,  setNames(as.list(c(if(use.rownames) {""},
                                                           names(data.frame(x,check.names = F)))),
                                                 MyFlexTable$header$dataset[1,]))
  MyFlexTable <- flextable::autofit(MyFlexTable)
  MyFlexTable <- flextable::hline_bottom(flextable::hline_top(flextable::border_remove(MyFlexTable),part='all',
                                                              border=officer::fp_border(color="black", width = 2)),
                                         border=officer::fp_border(color="black", width = 2),part='all')
  mydocFlex <- flextable::body_add_flextable(mydoc, MyFlexTable, align="center")
  if(!is.null(info)) mydocFlex <- officer::body_add_par(mydocFlex, info)
  print(mydocFlex,target= paste(file, ".docx", sep = ""))
}

#' Export a table to latex
#'
#' @description Exports a table to latex.
#' @param x A data.frame object
#' @param file Name of the file
#' @return Creates a .txt file with latex code for the table
#' @export
make_latex_table <- function(x, file){
  print(xtable::xtable(x), file=paste(file, ".txt", sep=""))
}

#' Export a table to excel
#'
#' @description Exports a table to Excel.
#' @param x A data.frame object
#' @param file Name of the file
#' @param info Footer for the table
#' @importFrom utils write.csv2
#' @return Creates a .csv file with the table
#' @export
make_csv_table <- function(x, file, info){
  write.csv2(data.frame(rbind(x, c(info, rep("", ncol(x)-1))), check.names=FALSE,
                        stringsAsFactors=FALSE), paste(file, ".csv", sep=""),
             row.names=TRUE)
}


#' Make a table from report
#'
#' @description Auxiliary function to create tables.
#' @param x A data.frame object
#' @param file Name of the file
#' @param type Type of file
#' @param info Footer for the table
#' @param ... Additional parameters passed to make_word_table
#' @return Creates a file with the table
#' @export
make_table<-function(x, file, type, info=NULL, ...){
  if(type=="csv") {make_csv_table(x, file, info)}
  if(type=="latex") {make_latex_table(x, file)}
  if(is.null(type) | type=="word") {make_word_table(x, file, info, ...)}
  message(paste0("Exported table as ", file))
}

#' Generic VarCorr function
#'
#' @description Extract Variance-Covariance Matrix.
#' @param x A model object
#' @param sigma Optional value used as a multiplier for the standard deviations
#' @param ... Further arguments passed to VarrCorr methods
#' @return A Variance-Covariance Matrix
#' @export
VarCorr<-function (x, sigma = 1, ...){
  UseMethod("VarCorr")
}

#' Set header names for word tables
#'
#' @description Internal function for make_word_table.
#' @param x A flextable object
#' @param args A names list with the header names
#' @importFrom utils tail
#' @return A flextable object with assigned header names
set_noms<-function (x, args){
  header_ <- x$header$dataset
  values <- as.list(tail(x$header$dataset, n = 1))
  args <- args[is.element(names(args), x$col_keys)]
  values[names(args)] <- args
  x$header$dataset <- rbind(header_[-nrow(header_), ], as.data.frame(values,
                                                                     stringsAsFactors = FALSE, check.names = FALSE))
  x
}

#' Generic function for creating reporting tables
#'
#' @description Generic function for creating reporting tables.
#' @param x An compatibleobject that can be summarized
#' @param ... further arguments passed to make_table
#' @return A data frame with the report table
#' @export
#' @examples
#' report(iris)  #Report of descriptive statistics
#' lm1 <- lm(Petal.Length ~ Sepal.Width + Species, data=iris)
#' report(lm1)   #Report of model
report<-function(x, ...){
  UseMethod("report")
}

#' Default function for report
#'
#' @description This is a default function for calling summary(x) on non-implemented classes.
#' @param x Any object without specific report function
#' @param ... further arguments passed to summary
#' @return A summary of the object
#' @export
report.default<-function(x, ...){
  warning(paste("Non-recognized class for report(). Returning summary(", deparse1(substitute(x)), ")", sep=""))
  return(summary(x, ...))
}

#' Report from numeric variable
#'
#' @description Creates a report table.
#' @param x A numeric variable
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @export
report.numeric<-function(x,...){
  report(data.frame(x), ...)
}

#' Report from categorical variable
#'
#' @description Creates a report table.
#' @param x A categorical variable
#' @param ... Further arguments passed to make_table
#' @return A data frame with the report table
#' @export
report.factor<-function(x,...){
  report(data.frame(x), ...)
}

#' Report tables of summary data
#'
#' @description Creates a report table ready for publication.
#' @param x A data.frame object
#' @param by Grouping variable for the report
#' @param remove.by Remove grouping variable from the report table?
#' @param file Name of the file to export the table
#' @param type Format of the file
#' @param digits Number of decimal places
#' @param digitscat Number of decimal places for categorical variables (if different to digits)
#' @param print Should the report table be printed on screen?
#' @param ... further arguments passed to make_table()
#' @return Returns a summary table of the data in publication-friendly format
#' @export
#' @examples
#' report(iris)
#' (reporTable<-report(iris, by="Species"))
#' class(reporTable)
report.data.frame <- function(x, by=NULL, remove.by=FALSE, file=NULL, type="word", digits=2, digitscat=digits, print=TRUE, ...){
  mtapply <- function(x, group, fun, simplify=TRUE){
    if(is.null(dim(x))) tapply(x, group, fun)
    else sapply(split(x, group, drop=TRUE), function(x) sapply(x, function(x) fun(x)), simplify = simplify)
  }
  if(!is.data.frame(x)){
    x <- data.frame(x)
  }
  if(is.null(by)) group <- factor(rep(1, nrow(x))) else group <- x[,by]
  if(remove.by) x <- x[,!names(x) %in% by, drop=FALSE]
  numeric_part <- x[, sapply(x, is.numeric), drop=FALSE]
  cat_part <- x[, sapply(x, is.factor), drop=FALSE]
  if(ncol(numeric_part)+ncol(cat_part) == 0) stop("No numeric or factor variables to describe. Please check that 'character' variables have been defined as 'factor'.")
  if(ncol(numeric_part)>0){
    #Numerical part
    num_res <- mtapply(numeric_part, group=droplevels(group), function(x) numeric_summary(x, digits=digits))
    format_num <- data.frame(Variable="", num_res)
    format_num[ ,1] <- unlist(lapply(names(numeric_part), function(x) c(paste(x, " (n)", sep=""), "  mean (sd)", "  median (q1, q3)")))
  } else format_num <- NULL
  if(ncol(cat_part)>0){
    #Categorical part
    cat_res <- matrix(unlist(mtapply(cat_part, group=droplevels(group), function(x) cat_summary(x, digits=digits), simplify=FALSE)), ncol=length(levels(interaction(group, drop=TRUE))))
    format_cat <- data.frame(Variable="", cat_res)
    format_cat[,1] <- unlist(lapply(names(cat_part), function(y) c(paste(y, " (n)", sep=""), paste("  ", levels(x[,y]), sep=""))))
    names(format_cat) <- make.names(c("Variable", levels(interaction(group, drop=TRUE))))
  } else format_cat <- NULL
  #Merge parts
  result <- rbind(format_num, format_cat)
  names(result) <- c(names(result)[1], paste(names(result)[-1], " (n = ", table(interaction(group, drop=TRUE)), ")", sep=""))
  if(!is.null(file)) make_table(result, file, type, use.rownames=FALSE)
  if(print){
    print(format(result, justify="left"), row.names=FALSE)
  } else{
    invisible(format(result, justify="left"))
  }
}

#' Auxiliary functions for report.data.frame
#' @description Internal function for numerical summary
#' @param x A vector of class numeric
#' @param digits Number of digits for rounding
#' @return Returns a summary of a numerical variable
numeric_summary <- function(x, digits){
  mean_sd <- function(x, digits){
    paste(round(mean(x, na.rm=TRUE), digits = digits), " (", round(sd(x, na.rm=TRUE), digits=digits), ")", sep="")
  }
  median_qq <- function(x, digits){
    paste(round(median(x, na.rm=TRUE), digits = digits), " (", paste(round(quantile(x, probs=c(0.25, 0.75), na.rm=TRUE), digits), collapse=", "), ")", sep="")
  }
  nsize <- function(x){
    sum(!is.na(x))
  }
  rbind(nsize(x), mean_sd(x, digits), median_qq(x, digits))
}

#' Auxiliary functions for report.data.frame
#' @description Internal function for categorical summary
#' @param x A vector of class character
#' @param digits Number of digits for rounding
#' @return Returns a summary of a categorical variable
cat_summary <- function(x, digits){
  c(sum(!is.na(x)), paste(table(x), " (", round(100*table(x)/sum(table(x)), digits), " %)", sep=""))
}

#' Auxiliary matrix paste function
#' @description Internal function for report.table
#' @param ... Matrices to paste
#' @param sep Separator for the paste function
#' @return Returns a matrix with the different matrices used as input pasted together
matrixPaste<-function (..., sep = rep(" ", length(list(...)) - 1)){
  theDots <- list(...)
  if (any(unlist(lapply(theDots, function(x) !is.character(x)))))
    stop("all matrices must be character")
  numRows <- unlist(lapply(theDots, nrow))
  numCols <- unlist(lapply(theDots, ncol))
  if (length(unique(numRows)) > 1 | length(unique(numCols)) >
      1)
    stop("all matrices must have the same dim")
  for (i in seq(along = theDots)) out <- if (i == 1)
    theDots[[i]]
  else paste(out, theDots[[i]], sep = sep[i - 1])
  matrix(out, nrow = numRows[1])
}

#' Plot of the coefficients of a model
#'
#' @description Creates a plot of the coefficients of a model.
#' @param coefs A vector with each coefficient
#' @param lwr.int A vector with the lower end of the CI
#' @param upper.int A vector with the upper end of the CI
#' @param offset Y-axis offset for the coefficients
#' @param coefnames Name for each variable
#' @param abline.pos Position for the vertical reference line
#' @param sorted Should the coefficients be sorted from highest to lowest?
#' @param reverse Should the order be reversed?
#' @param pch Type of point
#' @param xlim Limits of the X-axis
#' @param ylim Limits of the Y-axis
#' @param color Color for the points
#' @param ... Further arguments passed to axis()
#' @return A plot of the coefficients with their CI
#' @importFrom graphics abline lines plot.new plot.window points axis
#' @export
#' @examples
#' lm1 <- lm(Petal.Length ~ Sepal.Width + Species, data=iris)
#' a<-report(lm1)
#' oldpar <- par()
#' par(mar=c(4, 10, 3, 2))
#' #Coefplot calling plot.reportmodel
#' plot(a)
#' par(mar=oldpar$mar)  #Restore old margin values
#' #Manual coefplot
#' coefplot(coefs=c(1, 2, 3), lwr.int=c(0, 1, 2), upper.int=c(5, 6, 7), coefnames=c("A", "B", "C"))
coefplot <- function(coefs, lwr.int=coefs, upper.int=coefs, offset=0, coefnames=names(coefs), abline.pos=0, sorted=FALSE, reverse=FALSE, pch=16, xlim=c(min(lwr.int, na.rm=TRUE), max(upper.int, na.rm=TRUE)), ylim=c(1, length(coefs)), color="black", ...){
  color <- as.character(data.frame(color, coefs)[,1])
  if(is.null(coefnames)) coefnames <- 1:length(coefs)
  dat <- data.frame(coefs, lwr.int, upper.int, coefnames)
  if(sorted) dat <- dat[order(dat$coefs),]
  if(reverse) dat <- dat[(dim(dat)[1]):1,]
  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  axis(3, ...)
  axis(2, at=c(1:length(coefs)), las=1, labels=dat$coefnames, lwd=0)
  abline(v=abline.pos, lty=2)
  points(dat$coefs, (1:length(coefs))+offset, pch=pch, col=color)
  for(i in 1:length(coefs)){
    lines(x=c(dat$lwr.int[i], dat$upper.int[i]), y=c(i, i)+offset, col=color[i])
  }
}

#' Coefplot for reportmodel objects
#'
#' @description Creates a coefplot from the reportmodel object.
#' @param x A reportmodel object
#' @param ... Further arguments passed to coefplot
#' @return Returns a plot of each coefficient in the model with its 95% confidence interval
#' @export
#' @examples
#' lm1 <- lm(Petal.Length ~ Sepal.Width + Species, data=iris)
#' a<-report(lm1)
#' oldpar <- par()
#' par(mar=c(4, 10, 3, 2))
#' plot(a)   #Coefplot calling plot.reportmodel
#' par(mar=oldpar$mar)
plot.reportmodel<-function(x, ...){
  coefplot(x$coefficients, x$lwr.int, x$upper.int, ...)
}

#' K-fold cross-validation for models
#'
#' @description Returns cross-validated absolute fit values
#' @param formula An object of class "formula" describing the model to be validated
#' @param data A data frame containing the variables specified in formula argument
#' @param k Number of folds
#' @param fit_function Name of the model fitting function
#' @param metric Performance metric to estimate: RMSE, MSE, MAE  or AUC
#' @param predict.control Named list of arguments to pass to the predict function of the model
#' @param ... Further arguments passed to the model fitting function
#' @return Cross-validated values for the selected performance metric
#' @importFrom stats fitted predict
#' @export
#' @examples
#' cv_model(Petal.Length ~ Sepal.Width + Species, data=iris)
cv_model <- function(formula, data, k=5, fit_function="lm", metric=if(length(unique(data[,as.character(formula)[2]])) == 2) "AUC" else "RMSE", predict.control=list(NULL), ...){
  n <- nrow(data)
  fragmentos <- sample(1:k, n, replace=TRUE)
  sapply(1:k, function(x){
    mod_k <- get(fit_function)(formula, data=data[fragmentos != x,], ...)
    args <- append(list(object=mod_k, newdata = data[fragmentos == x,]), predict.control)
    if(fit_function == "brm"){
      prediccion_k <- do.call(fitted, args[!sapply(args, is.null)])
      get(metric)(prediccion_k[,1], data[,as.character(formula)[2]][fragmentos == x])
    } else{
      prediccion_k <- do.call(predict, args[!sapply(args, is.null)])
      get(metric)(prediccion_k, data[,as.character(formula)[2]][fragmentos == x])
    }
  })
}

#' Bootstrap optimism correction for models
#'
#' @description Returns optimism correction for absolute fit values
#' @param formula An object of class "formula" describing the model to be validated
#' @param data A data frame containing the variables specified in formula argument
#' @param B Number of bootstrap samples
#' @param fit_function Name of the model fitting function
#' @param metric Performance metric to estimate: RMSE, MSE, MAE  or AUC
#' @param predict.control Named list of arguments to pass to the predict function of the model
#' @param ... Further arguments passed to the model fitting function
#' @return Optimism correction values for the selected performance metric
#' @importFrom stats fitted predict
#' @export
#' @examples
#' boot_model(Petal.Length ~ Sepal.Width + Species, data=iris)
boot_model <- function(formula, data, B=200, fit_function="lm", metric=if(length(unique(data[,as.character(formula)[2]])) == 2) "AUC" else "RMSE", predict.control=list(NULL), ...){
  apparent_fit <- get(fit_function)(formula, data=data, ...)
  args_app <- append(list(object=apparent_fit, newdata = data), predict.control)
  prediction_app <- do.call(fitted, args_app[!sapply(args_app, is.null)])
  apparent_metric <- get(metric)(prediction_app, data[,as.character(formula)[2]])
  boot_samples <- replicate(B, data[sample(1:nrow(data), replace=TRUE), ], simplify = FALSE)
  boot_results <- do.call(rbind, lapply(boot_samples, function(x){
    mod_k <- get(fit_function)(formula, data=x, ...)
    args_boot <- append(list(object=mod_k, newdata = x), predict.control)
    args_orig <- append(list(object=mod_k, newdata = data), predict.control)
    if(fit_function == "brm"){
      prediccion_k_boot <- do.call(fitted, args_boot[!sapply(args_boot, is.null)])
      metric_boot <- get(metric)(prediccion_k_boot[,1], x[,as.character(formula)[2]])
      prediccion_k_orig <- do.call(fitted, args_orig[!sapply(args_orig, is.null)])
      metric_orig <- get(metric)(prediccion_k_orig[,1], data[,as.character(formula)[2]])
    } else{
      prediccion_k_boot <- do.call(predict, args_boot[!sapply(args_boot, is.null)])
      metric_boot <- get(metric)(prediccion_k_boot, x[,as.character(formula)[2]])
      prediccion_k_orig <- do.call(predict, args_orig[!sapply(args_orig, is.null)])
      metric_orig <- get(metric)(prediccion_k_orig, data[,as.character(formula)[2]])
    }
    data.frame(metric_boot=metric_boot, metric_orig=metric_orig)
  }))
  list(apparent_fit = apparent_metric, optimism=mean(boot_results$metric_boot - boot_results$metric_orig))
}


#' Mean squared error
#'
#' @description Estimates mean squared error from predicted and observed values
#' @param pred Numeric vector of predicted values
#' @param obs Numeric vector of observed values
#' @return Returns the MSE
#' @export
#' @examples
#' lm1 <- lm(Petal.Length ~ Sepal.Width + Species, data=iris)
#' pred1 <- fitted(lm1)
#' MSE(pred1, iris$Petal.Length)
MSE <- function(pred, obs){
  mean((pred-obs)^2)
}

#' Root mean squared error
#'
#' @description Estimates root mean squared error from predicted and observed values
#' @param pred Numeric vector of predicted values
#' @param obs Numeric vector of observed values
#' @return Returns the RMSE
#' @export
#' @examples
#' lm1 <- lm(Petal.Length ~ Sepal.Width + Species, data=iris)
#' pred1 <- fitted(lm1)
#' RMSE(pred1, iris$Petal.Length)
RMSE <- function(pred, obs){
  sqrt(MSE(pred, obs))
}

#' Mean absolute error
#'
#' @description Estimates mean absolute error from predicted and observed values
#' @param pred Numeric vector of predicted values
#' @param obs Numeric vector of observed values
#' @return Returns the MAE
#' @export
#' @examples
#' lm1 <- lm(Petal.Length ~ Sepal.Width + Species, data=iris)
#' pred1 <- fitted(lm1)
#' MAE(pred1, iris$Petal.Length)
MAE <- function(pred, obs){
  mean(abs(pred-obs))
}

#' ROC area under curve
#'
#' @description Estimates AUC from predicted and observed values
#' @param pred Numeric vector of predicted values
#' @param obs Numeric vector of observed values or factor with two levels
#' @return Returns the AUC
#' @export
#' @examples
#' glm1 <- glm(case ~ spontaneous + age, data=infert, family="binomial")
#' pred1 <- fitted(glm1)
#' AUC(pred1, infert$case)
AUC <- function(pred, obs){
  #Observed values
  ns <- table(obs)
  #Ranks in pred
  rank_pred <- rank(pred)
  #AUC
  max((sum(rank_pred[obs == names(ns)[2]]) - ns[2]*(ns[2]+1)/2)/(prod(ns)),
  (sum(rank_pred[obs == names(ns)[1]]) - ns[1]*(ns[1]+1)/2)/(prod(ns)))
}
