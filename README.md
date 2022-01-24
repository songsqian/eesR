# Environmental and Ecological Statistics with R

Code and data associated with [Environmental and Ecological Statistics with R](https://www.google.com/books/edition/Environmental_and_Ecological_Statistics/o6aKDQAAQBAJ?hl=en&gbpv=0) (2nd edition), CRC Press, 2016.

Example chapter -- chapter 11

R Update:

When using R-4.1.2 and newer, the self starter function (the initial value module) in Section 6.1.3 needs to be updated (page 235). In the initial value function in a self starter function must include `...` in its argument list:
```r
fplModelInit <- function(mCall, LHS, data, ...){
......
}
```
