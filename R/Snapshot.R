
.Snapshot <- setRefClass("Snapshot",
    fields=list(
      .debug="function", 
      .auto_display="logical",
      ## ranges
      .range="GRanges",
      .orig.range="GRanges",
      .zin="logical", 
      .pright="logical",
      ## data
      .data="data.frame",
      .data_dirty="logical",
      .initial_functions="SnapshotFunctionList",
      .current_function="character",
      .using_initial_functions="logical",
      ## annotation track 
      annTrack="ANY",
      ## more-or-less public
      functions="SnapshotFunctionList",
      files="BamFileList",
      view="SpTrellis"))


.Snapshot$methods(
    .message = function(fmt, ...) 
    {
        message(sprintf("Snapshot: %s", sprintf(fmt, ...)))
    },
    .stop=function(fmt, ...) 
    {
        stop(sprintf("Snapshot: %s", sprintf(fmt, ...)))
    },
    .initial_range=function() 
    {
        h <- scanBamHeader(.self$files[[1]])[["targets"]]
        if (!length(h))
            .stop("header of file '%s' contains no targets",
                  .self$files[[1]])
        h <- h[1]
        GRanges(names(h), IRanges(1, h))
    },

    .update_range=function(lim) {
        if (lim[2] < lim[1])
            .stop("The end of range must be greater than the start of the range.") 
        if (lim[1] >= start(.self$.orig.range)) {
            start(.self$.range) <- lim[1]
           .self$.data_dirty <- TRUE
       } else
           .stop("Please make sure the range argument are defining the regions within the limits of the original range.")
        if (lim[2] <= end(.self$.orig.range)) {
            end(.self$.range) <- lim[2]
            .self$.data_dirty <- TRUE
        } else
             .stop("Please make sure the range argument are defining the regions within the limits of the original range.")
        invisible()
    },
                  
    .update_data=function() 
    {
        .debug("update_data .current_function='%s'",
               .self$.current_function)
        .self$.data <-
            reader(.self$functions[[.self$.current_function]])(.self)
        .self$.data_dirty <- FALSE
        .self$view <-
            viewer(.self$functions[[.self$.current_function]])(.self)

        .debug("update_data view limits %.0f-%.0f",
               .self$view$get.x.limits()[1],
               .self$view$get.x.limits()[2])
        .self
    },

    .get.active_region=function() {
        'get the start and end of the active region'
        c(start(.self$.range), end(.self$.range))
    },

    .is.initial_function=function()
    {
        'check if initial reader/viwer function is currently in used:TRUE/FALSE'
        'assign result to .using_initial_functions'
         .self$.using_initial_functions <-
             any(.self$.current_function %in%
                 names(.self$.initial_functions)[1:2])
    },
                  
    .check_currentFunction=function(currentFunction) 
    {
        if (missing(currentFunction))
            currentFunction <- .self$.current_function

        lms <- limits(.self$functions[[currentFunction]])
        wd <- width(.self$.range)
        if (wd <= lms[1])
            .stop("image width (%.0f) < function limit (%.0f bps)",
                  wd, lms[1])
        ## FIXME: suggest to use togglefun to change function
        else if (wd > lms[2])
            .stop("image width (%.0f) > function limit (%.0f bps)",
                  wd, lms[2])
        invisible()
    },
                  
    .change_current_function=function(currentFunction) 
    {
        'Determine whether currentFunction should be change according to the
         size of the active region. This function is used by togglefun()' 
        'If yes, change .current_function and make .data_dirty TRUE'
        lms <- limits(.self$functions[[currentFunction]])
        wd <- .self$view$get.x.limits()[2] - .self$view$get.x.limits()[1]
        if (wd <= lms[1])
            .stop("image width (%.0f) < function limit (%.0f bps)",
                  wd, lms[1])
        ## FIXME: suggest to use togglefun to change function
        else if (wd > lms[2])
            .stop("image width (%.0f) > function limit (%.0f bps)",
                  wd, lms[2])
        .self$.current_function=currentFunction
        .self$.data_dirty <- TRUE
        invisible()
    },

    .zoom_in_xlim=function(){
        'get x limits for zoom in'
        lim <- .self$view$get.x.limits()
        center <- mean(lim)
        width <- (lim[2] - lim[1])/2
        if (width > 50)
            xlim <- c(max(start(.self$.orig.range), center - width/2),
                      min(end(.self$.orig.range), center + width/2))
        else xlim <- lim

        xlim[1] <- max(start(.self$orig.range), xlim[2])
    },
                  
    .zoom_out_xlim=function() {
        'get x limits for zoom out'
        lim <- .self$view$get.x.limits()
        center <- mean(lim)
        width <- diff(lim)
        xlim <-  c(max(start(.self$.orig.range), center-width),
                   min(end(.self$.orig.range), center+width))
    },
                  
    .pleft_xlim=function() {
        'get x limits for pan left'
        margin <- 50
        lim <- .self$view$get.x.limits()
        by <- 0.8 * diff(lim)
        xlim <- c(max(lim[1] - by, start(.self$.orig.range)),
                  max(lim[2] - by, start(.self$.orig.range) + margin))
        ## if xlim is between the gap of the limits of .self$range
        ## that of the trellis object limits (.self$view$trellis$orig.x.limits
        xlim <- c(min(end(.self$.orig.range)-margin, xlim[1]),
                  min(end(.self$.orig.range), xlim[2]))
    },
                  
    .pright_xlim=function() {
        'get x limits for pan right'
        margin <- 50
        lim <- .self$view$get.x.limits()
        by <- 0.8 * diff(lim)
        xlim <- c(min(lim[1]+by, end(.self$.orig.range) - margin),
                  min(lim[2]+by, end(.self$.orig.range)))
        ## if xlim is between the gap of the limits of .self$range
        ## that of the trellis object limits (.self$view$trellis$orig.x.limits
        xlim <- c(max(start(.self$.orig.range), xlim[1]),
                  max(start(.self$.orig.range)+margin, xlim[2]))
    },
                  
    .reset_active_range=function(xlim) {
       'determine wether to reset active range. used by pan and zoom out'
       win <- .self$view$trellis$orig.x.limits
       f1 <- xlim[1] < min(start(.self$.range), win[1])
       f2 <- xlim[2] > end(.self$.range, win[2])
       any(f1,f2)
    },
    
    .switch_ini_currentFunction=function(xlim) {
        'determine wether to switch viewer functions (TRUE/FALSE).
         used only when the current function is one of default functions
         (fine_coverage or coarse_coverage)'
        
        sw <- FALSE
        
        win <-(xlim[2] - xlim[1]) <
                   limits(.self$.initial_functions[["fine_coverage"]])[2]
        fine <- .self$.current_function == "fine_coverage"
        if (win) {
            # limits within fine_coverage limit and viewer is coarse
            if (!fine) sw <- TRUE
        } else {
            # limits over fine_coverage limit and viewer is fine
            if (fine) sw <- TRUE
        }
        
        return(sw)  
    },
                  
    .initialize_currentFunction=function() 
    {
        if (width(.self$.range) <
            limits(.self$.initial_functions[["fine_coverage"]])[2])
            currentFunction <- "fine_coverage"
        else currentFunction <- "coarse_coverage"
    },
                  
    initialize=function(..., functions=SnapshotFunctionList(), currentFunction,
                        .range, .auto_display=TRUE, .debug=FALSE, annTrack)
    {
        callSuper(...)
        .self$.debug <- if (.debug) .message else function(...) {}

        .self$.zin <- TRUE
        .self$.pright <- TRUE
        .self$.auto_display <- .auto_display
        tryCatch({
            for (f in as.list(.self$files))
                if (!isOpen(f)) open(f)
        }, error=function(err) {
            .stop("open BamFile failed: %s", conditionMessage(err))
        })

        .self$annTrack <-
            if (missing(annTrack)) NULL
            else annTrack
        
        .self$.range <- 
            if (missing(.range)) .initial_range()
            else .range
 
        .self$.orig.range <- .self$.range
        
        .self$.initial_functions <-
            SnapshotFunctionList(fine_coverage=.fine_coverage,
                                 coarse_coverage=.coarse_coverage,
                                 multifine_coverage=.multifine_coverage,
                                 multicoarse_coverage=.multicoarse_coverage)

        .self$functions <- c(.self$.initial_functions, functions)
         
            ## initialize current function
         if (!missing(currentFunction)) {
             if (!currentFunction %in% names(.self$functions))
                 .stop("'%s' not in SnapshotFunctionList",
                          currentFunction)
                 .self$.check_currentFunction(currentFunction)
         } else { 
                currentFunction <- .self$.initialize_currentFunction()
         }
        
        .self$.current_function <- currentFunction
        .self$.is.initial_function() # assign .self$using.initial_function
        .self$.data_dirty <- TRUE
        .self$.update_data()
        .self$display()
        .self
    },
                  
    set_range=function(range)
    {   'resetting the active range, called when setting zoom(..., range=)'
        'also used for determine the best fit SnapshotFunctions if the initial
         functions are in used.'
        # seqlevel must be the same
        if (!all(seqlevels(range) %in% seqlevels(.self$.range)))
           .stop("The seqlevel '%s' does not match that of the active data",
                 seqlevels(range))
        
        .self$.update_range(c(start(range), end(range)))
        .self$.is.initial_function()
        ## find appropriate reader/viewer if initial functions are in used
        if (.self$.using_initial_functions)
            .self$.current_function <- .self$.initialize_currentFunction()
        
        .self$.data_dirty <- TRUE
        .self$.update_data()
    },
                  
    display=function()
    {
        .debug("display")
        if (.data_dirty)
            .self$.update_data()
        print(.self$view$view())
    },

    toggle=function(zoom=FALSE, pan=FALSE, currentFunction)
    {
          .self$.debug("toggle: zoom %s; pan %s; fun %s",
                       if (.self$.zin) "in" else "out",
                       if (.self$.pright) "right" else "left",
                       .self$.current_function)
          if (zoom)
              .self$.zin <- !.self$.zin
          if (pan)
              .self$.pright <- !.self$.pright
          
          if (!missing(currentFunction)) {
              if (!currentFunction %in% names(.self$functions))
                  .stop("toggle unknown function '%s'", currentFunction)

              if (currentFunction != .self$.current_function) {
                 .self$.change_current_function(currentFunction)
                 if (.self$.data_dirty) {
                     lim <- .self$view$get.x.limits()
                     .update_range(lim)
                     .self$.update_data()
                 }
              }
          }
          .self
    },

    zoom=function()
    {
          .debug("zoom: %s", if (.self$.zin) "in" else "out")
          if (.self$.zin) {
              ## zoom in
              .self$.is.initial_function()
              if (.self$.using_initial_functions) {
                  # check if need to switch viewer
                  xlim <- .self$.zoom_in_xlim()
                  if (.self$.switch_ini_currentFunction(xlim)) {
                      range <- .self$.range
                      start(range) <- xlim[1]
                      end(range) <- xlim[2]
                      .self$set_range(range)
                  } else # if don't need to swith viewer
                      .self$view$zi()
              } else # if not using fine_coverage or coarse_coverage
              .self$view$zi()
          }
          else { ## zoom out
              xlim <- .self$.zoom_out_xlim()
              if (.reset_active_range(xlim)) { 
                  ## expend the active range and .update_data()
                  range <- .self$.range
                  start(range) <- xlim[1]
                  end(range) <- xlim[2]
                  #find appropriate read/viwer funcs
                  .self$set_range(range)
              }
              else
                  .self$view$zo()
          }
          #.self$.update_range()
          .self
    },                  

    pan=function() {
          .debug("pan: %s", if (.self$.pright) "right" else "left")
          if (.self$.pright) { ## shift right
              xlim <- .self$.pright_xlim()
              if (.reset_active_range(xlim)) {
                  .update_range(xlim)
                  .self$.update_data()
              }  
              else .self$view$right()
          }
          else { ## shift left
              xlim <- .self$.pleft_xlim()
              if (.reset_active_range(xlim)) {
                  .update_range(xlim)
                  .self$.update_data()
              }
              else .self$view$left()
          }    
          .self
    },
                  
    restore=function()
   {
       f1 <- start(.self$.range)==start(.self$.orig.range)
       f2 <- end(.self$.range)==end(.self$.orig.range)
       if (all(f1, f2))#original range is the same as active range
           .self$view$restore()
       else
          .self$set_range(.self$.orig.range) 
   }
)

## Constructors

setGeneric("Snapshot",
           function(files, range, ...) standardGeneric("Snapshot"))

setMethod(Snapshot, c("character", "GRanges"),
    function(files, range, ...)
{
    if (is.null(names(files)))
        names(files) <- basename(files)
    files <- BamFileList(files)
    .Snapshot$new(files=files, .range=range, ...)
})

setMethod(Snapshot, c("BamFileList", "GRanges"),
    function(files, range, ...)
{
    if (is.null(names(files))) 
        names(files) <- basename(sapply(files@listData, function(fl) path(fl)))

    .Snapshot$new(files=files, .range=range, ...)
})

setMethod(Snapshot, c("character", "missing"),
    function(files, range, ...)
{
    if (is.null(names(files)))
        names(files) <- basename(files)
    files <- BamFileList(files)
    .Snapshot$new(files=files, ...)
})

## accessors
if (is.null(getGeneric("files")))
    setGeneric("files", function(x, ...) standardGeneric("files"))
setMethod("files", "Snapshot", function(x)  x$files)


setGeneric("vrange", function(x, ...) standardGeneric("vrange"))
setMethod("vrange", "Snapshot", function(x) x$.range )

setGeneric("functions", function(x, ...) standardGeneric("functions"))
setMethod("functions", "Snapshot", function(x) x$functions)

setGeneric("annTrack", function(x, ...) standardGeneric("annTrack"))
setMethod("annTrack", "Snapshot", function(x) x$annTrack)

.getData <- function(x) x$.data
.currentFunction <- function(x) x$.current_function

if (is.null(getGeneric("view")))
  setGeneric("view", function(x, ...) standardGeneric("view"))

setMethod("view", "Snapshot", function(x) x$view)

## interface
setGeneric("togglez", function(x, ...) standardGeneric("togglez"))
setMethod("togglez", "Snapshot", function(x)
{
    x$toggle(zoom=TRUE)
    invisible(x)
})

setGeneric("togglep", function(x, ...) standardGeneric("togglep"))
setMethod("togglep", "Snapshot",  function(x)
{
    x$toggle(pan=TRUE)
    invisible(x)
})

setGeneric("togglefun", function(x, name, ...) standardGeneric("togglefun"))
setMethod("togglefun", "Snapshot",  function(x, name)
{
    if (!missing(name)) {
        x$toggle(currentFunction=name)
        invisible(x)
    }
})

setGeneric("zoom", function(x, range, ...) standardGeneric("zoom"))
setMethod("zoom", "Snapshot",  function(x, range)
{
    if (!missing(range))
        ## FIXME: must be able to tell whether .currentFunction is appropriate
        x$set_range(range) 
    else
        x$zoom()
    x$display()
    ## FIXME: invisible return TRUE on success, FALSE otherwise
})

setGeneric("pan", function(x, ...) standardGeneric("pan"))
setMethod("pan", "Snapshot", function(x)
{
    x$pan()
    x$display()
    ## FIXME: return TRUE on success, FALSE otherwise
})

## show
setMethod(show, "Snapshot", function(object) 
{
    cat("class:", class(object), "\n")
    with(object, {
        cat("file(s):", names(files), "\n")
        cat("Orginal range:",
            sprintf("%s:%d-%d", seqleveds(.orig.range), start(.orig.range),
                    end(.orig.range)), "\n")
        cat("active range:",
            sprintf("%s:%d-%d", seqlevels(.range), start(.range),
                    end(.range)), "\n")
        cat("zoom (togglez() to change):",
            if (.zin) "in" else "out", "\n")
        cat("pan (togglep() to change):",
            if (.pright) "right" else "left", "\n")
        cat("fun (togglefun() to change):",
            .current_function, "\n")
        cat(sprintf("functions: %s\n",
            paste(names(functions), collapse=" ")))
    })
    if (object$.auto_display)
        object$display()
})
