# Interactive plots with plotly

Every `ggvariant` plot function accepts `interactive = TRUE`, which
routes the plot through [plotly](https://plotly.com/r/) instead of
returning a static `ggplot` object. This is useful for sharing a plot
with collaborators who don’t use R: a plotly widget can be saved as a
standalone HTML file and opened in any browser.

`plotly` is a Suggests dependency, not a hard dependency of `ggvariant`
– the chunks on this page only run when `plotly` is installed on the
machine building the site.

``` r

library(ggvariant)

vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
variants <- read_vcf(vcf_file)
#> ℹ Reading VCF: example.vcf
#> ✔ Reading VCF: example.vcf [24ms]
#> 
#> Loaded 19 variant records across 7 chromosomes.
```

## Interactive lollipop plot

Hovering over a point shows the position, amino-acid change,
consequence, and sample.

``` r

plot_lollipop(variants, gene = "TP53", interactive = TRUE)
```

## Interactive oncoprint

``` r

plot_oncoprint(variants, top_n = 6, interactive = TRUE)
```

## Saving a widget for sharing

``` r

p <- plot_lollipop(variants, gene = "TP53", interactive = TRUE)
htmlwidgets::saveWidget(p, "TP53_lollipop.html")
```

[`htmlwidgets::saveWidget()`](https://rdrr.io/pkg/htmlwidgets/man/saveWidget.html)
is not a formal dependency of `ggvariant` – install it separately
(`install.packages("htmlwidgets")`) if you want to export a standalone
file this way.
