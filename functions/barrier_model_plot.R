##### BARRIER MODEL PLOTS
#https://haakonbakkagit.github.io/btopic128.html#62_Barrier_Model

plot.bmodel <- function(locs, spde, mesh, areapol, spmode = F, crs = NULL, range = 6, msd = 3, xl=NULL, yl=NULL,
                        plot.pol = T){
        
        library(fields)
        
        local.find.correlation = function(Q, location, mesh) {
                ## Vector of standard deviations
                sd = sqrt(diag(inla.qinv(Q)))
                
                ## Create a fake A matrix, to extract the closest mesh node index
                A.tmp = inla.spde.make.A(mesh=mesh, 
                                         loc = matrix(c(location[1],location[2]),1,2))
                
                ## Index of the closest node
                id.node = which.max(A.tmp[1, ])
                
                
                print(paste('The location used was c(', 
                            round(mesh$loc[id.node, 1], 4), ', ', 
                            round(mesh$loc[id.node, 2], 4), ')' ))
                
                ## Solve a matrix system to find the column of the covariance matrix
                Inode = rep(0, dim(Q)[1]) 
                Inode[id.node] = 1
                covar.column = solve(Q, Inode)
                # compute correaltions
                corr = drop(matrix(covar.column)) / (sd*sd[id.node])
                return(corr)
        }
        
        local.plot.field = function(field, mesh, xlim, ylim, ...){
                # Error when using the wrong mesh
                stopifnot(length(field) == mesh$n)
                
                # Choose plotting region to be the same as the study area polygon
                if (missing(xlim)) xlim = areapol@bbox[1, ] 
                if (missing(ylim)) ylim = areapol@bbox[2, ]
                
                # Project the mesh onto a 300x300 grid
                proj = inla.mesh.projector(mesh, xlim = xlim, 
                                           ylim = ylim, dims=c(300, 300))
                
                # Do the projection 
                field.proj = inla.mesh.project(proj, field)
                
                if (is.null(xl)) {
                        xl <- xlim
                }
                if(is.null(yl)) {
                        yl <- ylim
                }
                
                # Plot it
                image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
                           xlim = xl, ylim = yl, ...)  
        }
        
        
        if (spmode) {
                # spde1 NOT the barrier
                Q = inla.spde2.precision(spde, theta = c(log(range),log(msd)))

                corr = local.find.correlation(Q, loc = locs, mesh)

                local.plot.field(corr, mesh, zlim=c(0.1, 1))

                if (plot.pol) {
                        plot(areapol, add=T, col='grey', main = 'Stationary Model')
                }
                points(locs[1], locs[2], pch = 20, cex = 1, col = "yellow")
                
                title(main = 'Stationary Model')
        } else {
                Q = inla.rgeneric.q(b.model, "Q", theta = c(log(msd), log(range)))
                
                corr = local.find.correlation(Q, loc = locs, mesh)
                
                local.plot.field(corr, mesh, zlim=c(0.1, 1))
                if (plot.pol) {
                        plot(areapol, add=T, col='grey', main = 'Barrier Model')
                }
                points(locs[1], locs[2], pch = 20, cex = 1, col = "yellow")
                
                title(main = 'Barrier Model')
        }
        
}


# spde <- inla.spde2.pcmatern(
#         # Mesh and smoothness parameter
#         mesh = mesh, alpha = 2,
#         # P(practic.range < 0.3) = 0.5
#         prior.range = c(0.5, 0.01),
#         # P(sigma > 1) = 0.01
#         prior.sigma = c(1, 0.01)) 
# 
# pdf("teste.pdf")
# par(mfrow = c(2,2))
# s <- sample(1:nrow(sp.data$lyva), 30)
# for (i in 1:30) {
#         plot.bmodel(sp.data$lyva[s[i],1:2], mesh = mesh, spde = b.model, areapol = poly.barrier)
#         plot.bmodel(sp.data$lyva[s[i],1:2], mesh = mesh, spde = spde, areapol = poly.barrier, spmode = T)
# }
# dev.off()