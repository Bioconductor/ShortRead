setOldClass("trellis")

.SpTrellis <- setRefClass("SpTrellis",
    fields=list(
        trellis="trellis", 
        .debug_enabled="logical"))

.SpTrellis$methods(
    initialize=function(...) 
    {
        'initialize SpTrellis'
        callSuper(...)
        if (!is.null(.self$trellis))
            .self$ini.orig.x.limits(.self$get.x.limits())
        .self
    },
    get.trellis = function() 
    {
        'get trellis object'
        .self$trellis
    },
    view = function(window=NULL) 
    {
        'view trellis object (in a designated window)'
        margin=50
        if (is.null(window))
            .self$trellis
        else {
            if (window[2] < window[1]) {
                message("Invalid window")
                .self$trellis
            }
            else {
                if ((window[2] - window[1]) < margin)
                    window[2] <- window[1] + margin

                xlim <- c(max(window[1], .self$trellis$orig.x.limits[1]),
                          min(window[2], .self$trellis$orig.x.limits[2]))
                .self$trellis$x.limits <- xlim
                .self$trellis
            }
        }
    },
    get.x.limits = function() 
    {
        'get x.limits of a trellis object'
        .self$trellis$x.limits
    },
    get.y.limits = function() 
    {
        'get y.limits of a trellis object'
        .self$trellis$y.limits
    },
    ini.orig.x.limits = function(xlim) 
    {
        'set orig.x.limigs of the trellis object'
        .self$trellis$orig.x.limits <- xlim
        invisible(xlim)
    },
    set.x.limits = function(xlim) 
    {
        'set trellis object x.limits'
        if (xlim[1] < xlim[2])
            stop("Invalid x-axis limits")
        else .self$trellis$x.limits <- xlim
        invisible(xlim)
    },
    set.y.limits = function(ylim) 
    {
        'set trellis object .limits'
        if (ylim[1] < ylim[2])
            stop("Invalid x-axis limits")
        else .self$trellis$y.limits <- ylim
        invisible(ylim)
    },
    .debug = function(fmt, ...) {
        if (.self$.debug_enabled)
            message(sprintf("'.SpTrellis' %s", sprintf(fmt, ...)))
    },
        restore = function() {
        trellis$x.limits <<- trellis$orig.x.limits
        print(trellis)
    },
    zoomout = function() 
    {
        'zoom out'
        center <- mean(.self$trellis$x.limits)
        width <- (.self$trellis$x.limits[2]-.self$trellis$x.limits[1])
                 xlim <- c(center - width, center + width)
         
        xlim[1] <- max(.self$trellis$orig.x.limits[1], xlim[1])
        xlim[2] <- min(.self$trellis$orig.x.limits[2], xlim[2])
        .self$trellis$x.limits <- xlim
        .debug("current x limits [%.0f, %.0f]", .self$get.x.limits()[1],
                        .self$get.x.limits()[2])
    },
    zoomin = function() {
        'zoom in 50%, change x.limits'
         center <- mean(.self$trellis$x.limits)
         width <- (.self$trellis$x.limits[2] - .self$trellis$x.limits[1])/2
         if (width > 1)
             xlim <- c(center - width/2, center + width/2)
         else xlim <- .self$trellis$x.limits
         .self$trellis$x.limits <- xlim
         .debug("current x limits [%.0f, %.0f]", 
                .self$get.x.limits()[1], .self$get.x.limits()[2])
    },
    shiftl = function(by=NULL) 
    {
        'shift x.limits 80% to the left'
        margin <- 50
         
       if (is.null(by)) ## 80% to the left
           by <- 0.8 * (.self$trellis$x.limits[2] - .self$trellis$x.limits[1])
           .debug("shift left for %.0f bps", by)
           .self$trellis$x.limits[1] <- max(.self$trellis$x.limits[1] - by,
                                       .self$trellis$orig.x.limits[1])
           .self$trellis$x.limits[2] <- max(.self$trellis$x.limits[2] - by,
                                        .self$trellis$orig.x.limits[1] + margin)
           .debug("current x limits [%.0f, %.0f]", .self$get.x.limits()[1],
                  .self$get.x.limits()[2])
    },
    shiftr = function(by=NULL) 
    {
        'shift x.limits 80% to the right'
         margin <- 50 #pbs
         if (is.null(by)) ## 80% to the left
             by <- 0.8 * (.self$trellis$x.limits[2] - .self$trellis$x.limits[1])
         .debug("shift right for %s bps", by)
         .self$trellis$x.limits[1] <- min(.self$trellis$x.limits[1] + by,
                                    .self$trellis$orig.x.limits[2] - margin)
         .self$trellis$x.limits[2] <- min(.self$trellis$x.limits[2] + by,
                                     .self$trellis$orig.x.limits[2])
         .debug("current x limits [%.0f, %.0f]", .self$get.x.limits()[1],
                .self$get.x.limits()[2])
    },
    current = function() 
    {
        print(.self$trellis)
    },
    restore = function()
    {
        .self$trellis$x.limits <- .self$trellis$orig.x.limits
        ## .self$trellis
    }

                   )

SpTrellis <- function(trellis, debug_enabled=FALSE) 
{   
    .SpTrellis$new(trellis = trellis, .debug_enabled=debug_enabled)

}
            
                    
                   
                   
