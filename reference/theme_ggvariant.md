# ggvariant ggplot2 theme

A clean, publication-ready theme based on `theme_minimal`. Applied
automatically by all `ggvariant` plot functions; export it to customise
further.

## Usage

``` r
theme_ggvariant(base_size = 12, base_family = "")
```

## Arguments

- base_size:

  Base font size in pt. Default `12`.

- base_family:

  Base font family. Default `""` (system sans-serif).

## Value

A `ggplot2` theme object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_ggvariant()

```
