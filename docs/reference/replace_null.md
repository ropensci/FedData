# Replace NULLs

Replace all the empty values in a list

## Usage

``` r
replace_null(x)
```

## Arguments

- x:

  A list

## Value

A list with NULLs replaced by NA

## Examples

``` r
list(a = NULL, b = 1, c = list(foo = NULL, bar = NULL)) %>% replace_null()
#> $a
#> [1] NA
#> 
#> $b
#> [1] 1
#> 
#> $c
#> $c$foo
#> NULL
#> 
#> $c$bar
#> NULL
#> 
#> 
```
