# Remove duplicate points, keeping just one per cell
# This is to be used with Presence-Absence data
# data = SpatialPointsDataFrame
# env = base layer to be used
dup.cells <- function(data, env) {
  
  #Create a clean raster
  r <- env
  r[] <- NA
  
  #Put 1 in presence cells
  for (i in 1:nrow(data)) {
    r[terra::cellFromXY(r, matrix(data@coords[i,], ncol = 2))] <- unlist(ifelse(
      is.na(r[terra::cellFromXY(r, matrix(data@coords[i,], ncol = 2))]),
      data$presence[i],
      r[terra::cellFromXY(r, matrix(data@coords[i,], ncol = 2))] + data$presence[i]))
  }
  
  r[r > 1] <- 1
  
  #Convert raster to dataframe
  data <- as.points(r, values=TRUE, na.rm=TRUE, na.all=FALSE)
  data <- as.data.frame(data, geom = "XY")
  
  data <- SpatialPointsDataFrame(
    data[,2:3], data.frame(presence = data[,1]),
    proj4string = CRS(crs(env, proj = T))
  )
  
  data
  
}

# Functions for easier ploting
plot.res <- function(var = NULL, # Which variable to plot
                     obj = pred.imodel, # The object holding predictions
                     metric = "mean", # Which metric (mean, sd, etc.)
                     po = F, ab = F, pa = F, # Plot points?
                     xlim = NULL, ylim = NULL){ # To zoom in a region
  
  sca <- function(...){
    scale_fill_gradientn(
      colours = rev(c("#A84C00", "#D97D27", "#F5BD44", "#FFD561", "#FFF291",
                      "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", 
                      "#313695")),
      limits = range(...),
      #labels = c("Low", "", "", "", "High"),
      na.value = "#572701", name = "Value")
  }
  
  if (metric != "iqr") {
    if (is.null(var)) {
      p <- ggplot()+gg(obj, aes(fill = eval(parse(text = metric))))+
        sca(obj[[metric]])+
        coord_equal(ylim = ylim, xlim = xlim)
    } else{
      p <- ggplot()+gg(obj[[var]], aes(fill = eval(parse(text = metric))))+
        sca(obj[[var]][[metric]])+
        coord_equal(ylim = ylim, xlim = xlim)
    }
    if (po) {p <- p+gg(po.pts, color = "red", size = .5)}
    if (ab) {p <- p+gg(ab.pts, color = "blue", aes(size = abundance), alpha = .1)}
    if (pa) {p <- p+gg(pa.pts, aes(color = as.factor(presence)))}
  } else{
   
    
    if (is.null(var)) {
      vals <- c(obj[["q0.025"]], obj[["mean"]], obj[["q0.975"]])
    } else {
      vals <- c(obj[[var]][["q0.025"]], obj[[var]][["mean"]], obj[[var]][["q0.975"]])
    }
    
    n.obj <- as.data.frame(obj[[var]][c("q0.025", "mean", "q0.975")])
    n.obj <- tidyr::pivot_longer(n.obj, cols = 1:3, names_to = "metric",
                                 values_to = "value")
    
    p <- ggplot(n.obj)+
      geom_tile(aes(x = x, y = y, fill = value))+
      sca(vals)+
      coord_equal(ylim = ylim, xlim = xlim) +
      facet_wrap(~metric)
    
    if (po) {p <- p+gg(po.pts, color = "red", size = .5)}
    if (ab) {p <- p+gg(ab.pts, color = "blue", aes(size = abundance), alpha = .1)}
    if (pa) {p <- p+gg(pa.pts, aes(color = as.factor(presence)))}
    
  }
  p
}

# Plot all covariates or a set of predictions
plot.subset <- function(var = NULL, # vector of predictions to plot,
                        # if NULL, covariates will be ploted
                        obj = pred.imodel, # The object holding predictions
                        metric = "mean", # Which metric (mean, sd, etc.)
                        po = F, ab = F, pa = F, # Plot points?
                        xlim = NULL, ylim = NULL){ # To zoom in a region
  
  library(patchwork)
  
  if (is.null(var)) {
    var <- obj$names.fixed
    var <- var[-grep("tercept", var)]
  }
  
  sca <- function(...){
    scale_fill_gradientn(
      colours = rev(c("#A84C00", "#D97D27", "#F5BD44", "#FFD561", "#FFF291",
                      "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", 
                      "#313695")),
      limits = range(...),
      #labels = c("Low", "", "", "", "High"),
      na.value = "#572701", name = "Value")
  }
  
  plot.list <- list()
  
  for (i in 1:length(var)) {
    p <- ggplot()+gg(obj[[var[i]]], aes(fill = eval(parse(text = metric))))+
      sca(obj[[var[i]]][[metric]])+
      coord_equal(ylim = ylim, xlim = xlim)
    if (po) {p <- p+gg(po.pts, color = "red", size = .5)}
    if (ab) {p <- p+gg(ab.pts, color = "blue", aes(size = abundance), alpha = .1)}
    if (pa) {p <- p+gg(pa.pts, aes(color = as.factor(presence)))}
    plot.list[[i]] <- p
  }
  
  for (i in 1:length(var)) {
    if (i == 1) {
      fp <- plot.list[[1]]
    } else{
      fp <- fp + plot.list[[i]]
    }
  }
  
  fp
}

# Plot temperature SPDE
plot.temp <- function(x, knots = TRUE, tr = NULL, mode = "d1"){
  dat <- data.frame(
    x = round(x$summary.random[[1]]$ID, 2),
    mean = x$summary.random[[1]]$mean,
    up = x$summary.random[[1]]$`0.975quant`,
    lo = x$summary.random[[1]]$`0.025quant`
  )
  
  if (mode == "d1") {
    if (knots) {
      dat$x <- round(d1mesh.st$mid, 2)
    }
  }
  
  p <- ggplot(dat) +
    geom_line(aes(x = x, y = mean)) +
    geom_ribbon(aes(x = x, ymin = lo, ymax = up), alpha = .4) +
    scale_x_continuous(breaks = seq(min(dat$x), max(dat$x), by = 0.5)) +
    geom_hline(yintercept = 0, color = "grey50")+
    theme_classic()
  
  if (!is.null(tr)) {
    p <- p + geom_vline(xintercept = tr[1], color = "blue") +
      geom_vline(xintercept = tr[2], color = "blue")
  }
  
  p
}

# 
# df.teste <- data.frame(env.e = seq(-3.7, 1.2, by = 0.1)^2,
#                        env.e2 = seq(-3.7, 1.2, by = 0.1))
# 
# pred.teste <- predict(imodel, df.teste,
#                       ~ tempmax_sq_eval(env.e) + tempmax_eval(env.e2))
# 
# plot(pred.teste$mean ~ seq(-3.7, 1.2, by = 0.1))
# 
# pred.teste.glm <- predict(teste1,
#                           data.frame(
#                             tempmax = seq(-3.7, 1.2, by = 0.1),
#                             salinitymean = mean(data.teste$salinitymean, na.rm = T),
#                             chlomean = mean(data.teste$chlomean, na.rm = T),
#                             ph = mean(data.teste$ph, na.rm = T)
#                           ))

delta.metrics <- function(df, metric, all = F){
  
  if (all(all != FALSE)) {
    if (length(all) > 1) {
      if (!all(all %in% names(df))) {
        stop("Metrics should be one of: ", paste(names(df), collapse = " "))
      }else{
        metrics <- all
      }
      
    } else{
      metrics <- names(df)[-grep("model", names(df))]
    }
    
    for (i in metrics) {
      
      if (grepl("auc", i) | grepl("tss", i) | grepl("boyce", i)) {
        dec = TRUE
      } else{
        dec = FALSE
      }
      
      if (i == metrics[1]) {
        dfb <- df[,c("model", i)]
        
        dfb <- dfb[order(dfb[,2], decreasing = dec),]
        
        dfb$delta <- dfb[,2] - dfb[1,2]
        names(dfb) <- c("model", i, paste(i, "delta", sep = "_"))
      } else{
        dfc <- df[,c("model", i)]
        
        dfc <- dfc[order(dfc[,2], decreasing = dec),]
        
        dfc$delta <- dfc[,2] - dfc[1,2]
        names(dfc) <- c("model", i, paste(i, "delta", sep = "_"))
        
        dfb <- cbind(dfb, dfc)
      }
    }
    
    df <- dfb
    
  } else {
    if (metric %in% names(df) == 0) {
      stop("Metric should be one of: ", paste(names(df), collapse = " "))
    }
    
    if (grepl("auc", metric)) {
      dec = TRUE
    } else{
      dec = FALSE
    }
    
    df <- df[,c("model", metric)]
    
    df <- df[order(df[,2], decreasing = dec),]
    
    df$delta <- df[,2] - df[1,2]
  }
  
  return(df)
  
}


# Check if all data is valid in the extended dataset
check.extend <- function(data.po, data.pa, data.ips, lays){
  po <- all(!is.na(extract(lays, coordinates(data.po))[,1]))
  pa <- all(!is.na(extract(lays, coordinates(data.pa))[,1]))
  ip <- all(!is.na(extract(lays, coordinates(data.ips))[,1]))
  
  if (all(c(po, pa, ip))) {
    cat("All ok! \n")
    return(invisible(NULL))
  } else {
    warning("One of the data have problems! \n")
    return(c("PO", "PA", "IPS")[!c(po, pa, ip)])
  }
}
