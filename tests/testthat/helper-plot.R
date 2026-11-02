# Evaluates 'code' with a pdf device on a temporary file, so that the plots of the tests do not go
# to Rplots.pdf, and returns the value of 'code' with the size of the file as the attribute
# "pdf_size"
with_pdf <- function(code) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  res <- tryCatch(force(code), finally = grDevices::dev.off())
  size <- file.size(f)
  unlink(f)
  attr(res, "pdf_size") <- size
  res
}

# Design B of test-ms_power.R (non-binding futility bounds by beta spending)
plot_design_b <- function() {
  te <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1)
  ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)), test = te,
           events = 120)
}
