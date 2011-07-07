
.Snapshot <- setRefClass("Snapshot",
    fields=list(
      .debug="function", 
      .auto_display="logical",
      ## ranges
      .range="GRanges", 
      .zin="logical", 
      .pright="logical",
      ## data
      .data="data.frame",
      .data_dirty="logical",
      .initial_functions="SnapshotFunctionList",
      .current_function="character",
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
    .update_range=function() 
    {
        lim <- pmax(1L, .self$view$get.x.limits())
        if (lim[1] < start(.self$.range)) {
            start(.self$.range) <- lim[1]
            .self$.data_dirty <- TRUE
        }
        if (lim[2] > end(.self$.range)) {
            end(.self$.range) <- lim[2]
            .self$.data_dirty <- TRUE
        }
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
        'check whether currentFunction should be change according to the window'
        'if yes, change .current_function and make .data_dirty TRUE'
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

            .self$.initial_functions <-
                SnapshotFunctionList(fine_coverage=.fine_coverage,
                                     coarse_coverage=.coarse_coverage,
                                     multifine_coverage=.multifine_coverage)

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
            .self$.data_dirty <- TRUE
            .self$.update_data()
            .self$display()
            .self
    },
                  
    set_range=function(range)
    {
        .self$.range <- range
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
                     #.self$.update_range()
                     # just update .range here instead of using .update_range()
                     lim <- .self$view$get.x.limits()
                     start(.self$.range) <- lim[1]
                     end(.self$.range) <- lim[2]
                     .self$.update_data()
                 }
                 #.self$data_dirty <- TRUE
                 #.self$update_data()
                
              }
          }
          .self
    },

    zoom=function()
    {
          .debug("zoom: %s", if (.self$.zin) "in" else "out")
          if (.self$.zin)
              .self$view$zoomin()
          else
              .self$view$zoomout()
          .self$.update_range()
          .self
    },                  

    pan=function() {
          .debug("pan: %s", if (.self$.pright) "right" else "left")
          if (.self$.pright)
              .self$view$shiftr()
          else
              .self$view$shiftl()
          .self$.update_range()
          .self
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
