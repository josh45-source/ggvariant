## R CMD check results

0 errors | 0 warnings | 1 note

One NOTE was observed locally:

  * checking for future file timestamps ... NOTE
    unable to verify current time

This is an environmental NOTE from a failed outbound query to a time
verification service; it is not caused by anything in the package and did
not reproduce consistently across repeated local check runs.

## This is a patch release

ggvariant 0.1.1 fixes two documented-but-unimplemented arguments discovered
after the 0.1.0 release:

  * `plot_variant_spectrum()`'s `context` and `genome` arguments, and
    `read_vcf()`'s `info_fields` argument, were fully documented (including
    a claimed BSgenome/Biostrings dependency for `genome`) but never read by
    their function bodies. A call such as
    `plot_variant_spectrum(x, context = TRUE)` silently returned the default
    6-class plot with no error or warning, while the caller believed they
    had requested a 96-trinucleotide-context plot. Each now aborts with an
    informative error until the feature is implemented.

This release also removes several stray scratch files that had been
accidentally committed to the package tarball.
