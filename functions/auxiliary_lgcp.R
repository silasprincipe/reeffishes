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
    r[cellFromXY(r, data@coords[i,])] <- ifelse(
      is.na(r[cellFromXY(r, data@coords[i,])]),
      data$presence[i],
      r[cellFromXY(r, data@coords[i,])] + data$presence[i])
  }
  
  r[r > 1] <- 1
  
  #Convert raster to dataframe
  data <- rasterToPoints(r)
  
  data <- SpatialPointsDataFrame(
    data[,1:2], data.frame(presence = data[,3]),
    proj4string = crs(env)
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
    x = round(x$summary.random$tempmax$ID, 2),
    mean = x$summary.random$tempmax$mean,
    up = x$summary.random$tempmax$`0.975quant`,
    lo = x$summary.random$tempmax$`0.025quant`
  )
  
  if (mode == "d1") {
    if (knots) {
      dat$x <- d1mesh.st$mid
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
      
      if (grepl("auc", i)) {
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
