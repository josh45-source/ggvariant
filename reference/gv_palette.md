# ggvariant colour palettes

Access the built-in colour palettes used by `ggvariant` plot functions.

## Usage

``` r
gv_palette(type = c("consequence", "spectrum", "domain"), n = 8L)
```

## Arguments

- type:

  One of `"consequence"` (default), `"spectrum"`, or `"domain"`.

- n:

  Integer. For `"domain"`, the number of colours to generate.

## Value

A named character vector of hex colour codes.

## Details

The `"consequence"` palette covers the most common VEP/SnpEff/MAF
consequence terms; any other term is standardised to an `"Other"` bucket
(its own entry in this palette) by plot functions rather than being
dropped.

## Examples

``` r
gv_palette("consequence")
#>     missense_variant             missense    Missense_Mutation 
#>            "#FD8D3C"            "#FD8D3C"            "#FD8D3C" 
#>   synonymous_variant           synonymous               Silent 
#>            "#74C476"            "#74C476"            "#74C476" 
#>   frameshift_variant  frameshift_deletion frameshift_insertion 
#>            "#E31A1C"            "#E31A1C"            "#FC4E2A" 
#>          stop_gained    Nonsense_Mutation            stop_lost 
#>            "#800026"            "#800026"            "#BD0026" 
#>  splice_site_variant          Splice_Site    inframe_insertion 
#>            "#9E9AC8"            "#9E9AC8"            "#6BAED6" 
#>     inframe_deletion         In_Frame_Del         In_Frame_Ins 
#>            "#2171B5"            "#2171B5"            "#6BAED6" 
#>  5_prime_UTR_variant  3_prime_UTR_variant       intron_variant 
#>            "#BCBDDC"            "#DADAEB"            "#D9D9D9" 
#>                  SNV             deletion            insertion 
#>            "#BDBDBD"            "#41AB5D"            "#A1D99B" 
#>                  MNV                Other 
#>            "#F768A1"            "#7F7F7F" 
gv_palette("spectrum")
#>       C>A       C>G       C>T       T>A       T>C       T>G 
#> "#03BDEE" "#010101" "#E52926" "#CAC9C9" "#A0CE63" "#ECC7C5" 
```
