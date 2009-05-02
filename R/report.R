## PDF

.report_pdf_do <-
    function(src, dest, symbolValues) 
{
    if (!file.exists(dirname(dest)))
        .throw(SRError("Input/Output",
                       "'dest' directory '%s'\n  does not exist",
                       dirname(dest)))
    if (file.exists(dest))
        .throw(SRError("Input/Output",
                       "'dest' file '%s'\n  already exists",
                       dest))
    tmpdir <- tempfile()
    if (!dir.create(tmpdir))
        .throw(SRError("Input/Output",
                       "failed to create temporary directory '%s'",
                       tmpdir))
    cwd <- setwd(tmpdir)
    on.exit(setwd(cwd))

    tmpfile <- file.path(tmpdir, basename(src))
    copySubstitute(src, tmpfile, symbolValues)
    texFile <- Sweave(tmpfile)
    tools::texi2dvi(texFile, pdf=TRUE)
    o_pdfFile <- sub(".tex$", ".pdf", texFile)
    ok <- file.copy(o_pdfFile, dest)
    if (!ok)
        .throw(SRError("Input/Output",
                       "failed to copy '%s'\n  to '%s'",
                       o_pdfFile, dest))
    dest
}

setMethod(.report_pdf, "character",
          function(x, dest, type, ...)
{
    src <- system.file("template", "qa_solexa.Rnw",
                       package="ShortRead")
    if (.Platform$OS.type == "windows")
      x <- gsub("\\\\", .Platform$file.sep, x)
    symbolValues <- list(QA_SAVE_FILE=x)
    .report_pdf_do(src, dest, symbolValues)
})

## HTML

.report_html_do <-
    function(destDir, sections, values, ...)
{
    htmlFile <- file.path(destDir, "index.html")
    cssFile <- "QA.css"
    biocFile <- "bioclogo-small.jpg"
    values <-
        c(list(CSS=cssFile, DATE=date(),
               VERSION=packageDescription("ShortRead")$Version),
          values)
    toConn <- file(htmlFile, "w")
    for (sec in sections) {
        fromConn <- file(sec, open="r")
        copySubstitute(sec, toConn, values)
        close(fromConn)
    }
    close(toConn)
    file.copy(system.file("template", cssFile, package="ShortRead"),
              file.path(destDir, cssFile))
    file.copy(system.file("template", "image", biocFile,
                          package="ShortRead"),
              file.path(destDir, "image", biocFile))
    htmlFile
}

.html_NA <- function() "<pre>NA</pre>"

.html_img <-
    function(dir, file, fig, ...)
{
    jpegFile <- paste(file, "jpg", sep=".")
    pdfFile <- paste(file, "pdf", sep=".")
    imgDir <- file.path(dir, "image")
    if (!file.exists(imgDir))
        dir.create(imgDir)

    jpeg(file.path(imgDir, jpegFile), ...)
    print(fig)
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)
    print(fig)
    dev.off()

    hwriteImage(file.path(".", "image", jpegFile),
                link=file.path(".", "image", pdfFile))
}

.htmlReadQuality <-
    function(dir, file, qa, type="read", ...)
{
    df <- qa[["readQualityScore"]]
    .html_img(dir, file,
              .plotReadQuality(df[df$type==type,]),
              ...)
}

.htmlReadOccur <-
    function(dir, file, qa, type="read", ...)
{
    df <- qa[["sequenceDistribution"]]
    .html_img(dir, file,
              .plotReadOccurrences(df[df$type==type,], cex=.5),
              ...)
}
