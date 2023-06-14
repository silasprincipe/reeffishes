### Get response curves from inlabru predictions
get.resp.curves <- function(fit, vars, mode = NULL, samp = 100,
                            int = FALSE, used = TRUE){
  
  if (class(vars) == "formula") {
    vars <- all.vars(vars)
  }
  
  if (isTRUE(used)) {
    env.e <- rast(c("data/env/ready_layers/tempmean.tif",
                    "data/env/ready_layers/tempmax.tif",
                    "data/env/ready_layers/salinitymean.tif",
                    "data/env/ready_layers/chlomean.tif",
                    "data/env/ready_layers/silicatemax.tif",
                    "data/env/ready_layers/ph.tif",
                    "data/env/ready_layers/windspeed.tif",
                    "data/env/ready_layers/distcoast.tif"))
    bvals <- rbind(extract(env.e, coordinates(ips)),
                   extract(env.e, coordinates(po.pts)),
                   extract(env.e, coordinates(pa.pts)))
    bvals <- bvals[,vars]
    min.vals <- apply(bvals, 2, min)
    max.vals <- apply(bvals, 2, max)
    mea.vals <- apply(bvals, 2, mean)
  }else{
    env.data <- env[[vars]]
    min.vals <- global(env.data, "min", na.rm=T)[,1]
    max.vals <- global(env.data, "max", na.rm=T)[,1]
    mea.vals <- global(env.data, "mean", na.rm=T)[,1]
  }
  
  pred.list <- list()
  
  resp.data <- matrix(nrow = 100, ncol = length(vars))
  
  for (z in 1:length(vars)) {
    
    resp.data[,] <- rep(mea.vals, each = 100)
    
    cat("Predicting response curve for", vars[z], "\n")
    
    resp.data[,z] <- seq(min.vals[z], max.vals[z], length.out = 100)
    
    if (!is.null(mode)) {
      if (mode == "cloglog") {
        pred.resp <- predict(fit,
                             NULL,
                             ~ eval(parse(text =
                                            paste0("1-exp(-exp(",
                                                   paste0(vars,
                                                          "_eval(resp.data[,",
                                                          1:length(vars),
                                                          "])", collapse = "+"),
                                                   ifelse(int, "+intercept_pa_latent", ""),
                                                   "))")
                             )), n.samples = samp
        )
      } else{
        pred.resp <- predict(fit,
                             NULL,
                             ~ eval(parse(text =
                                            paste0(mode, "(",
                                                   paste0(vars,
                                                          "_eval(resp.data[,",
                                                          1:length(vars),
                                                          "])", collapse = "+"),
                                                   ifelse(int, "+Intercept_latent", ""),
                                                   ")")
                             )), n.samples = samp
        )
      }
    } else{
      pred.resp <- predict(fit,
                           NULL,
                           ~ eval(parse(text =
                                          paste0(mode, "(",
                                                 paste0(vars,
                                                        "_eval(resp.data[,",
                                                        1:length(vars),
                                                        "])", collapse = "+"),
                                                 ifelse(int, "+Intercept_latent", ""),
                                                 ")")
                           )), n.samples = samp
      )
    }
    
    pred.resp$base <- resp.data[,z]
    
    pred.list[[z]] <- pred.resp
    
  }
  
  if (!is.null(mode)) {
    tp <- mode
  } else{
    tp <- "lin"
  }
  
  names(pred.list) <- vars
  
  return(structure(pred.list, class = c("respCur", tp, class(pred.list))))
}

# Old version, multiple plots
# plot.respCur <- function(x) {
#   
#   plot.list <- list()
#   
#   sc <- ifelse("lin" %in% class(x), "Linear predictor",
#                ifelse("exp" %in% class(x), "Relative Occurrence Rate", "Probability"))
#   
#   for (z in 1:length(x)) {
#     
#     df <- x[[z]]
#     
#     plot.list[[z]] <- ggplot(df) +
#       geom_line(aes(x = base, y = mean)) +
#       geom_ribbon(aes(x = base, ymin = q0.025, ymax = q0.975), alpha = .4, fill = "#52638E") +
#       scale_x_continuous(expand = c(0,0))+
#       theme_classic() + ylab(sc) + xlab(names(x)[z])
#     
#   }
#   
#   print(eval(parse(text = paste("plot.list[[", 1:length(plot.list), "]]", collapse = "+"))))
#   
#   return(invisible(NULL))
#   
# }

plot.respCur <- function(x, free = NULL, mode = "both") {
  
  sc <- ifelse("lin" %in% class(x), "Linear predictor",
               ifelse("exp" %in% class(x), "Relative Occurrence Rate", "Probability"))
  
  if (is.null(free)) {
    if (sc == "Probability") {
      free = F
    } else {
      free = T
    }
  }
  
  x <-  lapply(names(x), function(z){cbind(x[[z]], var = z)})
  
  dat <- do.call("rbind", x)
  
  if (mode != "both") {
    if (mode == "median") {
      p <- ggplot(dat) +
        geom_line(aes(x = base, y = q0.5))
    } else {
      p <- ggplot(dat) +
        geom_line(aes(x = base, y = mean))
    }
  } else {
    p <- ggplot(dat) +
      geom_line(aes(x = base, y = mean), color = "#233461", linetype = "dashed") +
      geom_line(aes(x = base, y = q0.5))
  }
  
  p <- p + 
    geom_ribbon(aes(x = base, ymin = q0.025, ymax = q0.975), alpha = .4, fill = "#52638E") +
    scale_x_continuous(expand = c(0,0))+
    theme_classic() + ylab(sc) + xlab(NULL) +
    theme(strip.background = element_blank()) +
    facet_wrap(~var, scales = ifelse(free, "free", "free_x"))
    
  
  print(p)
  
  return(invisible(NULL))
  
}
