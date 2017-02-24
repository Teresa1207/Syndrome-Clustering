#setwd("~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ServerResults")
#source("Results_Serv.R")
require(xtable)
results = data.frame(matrix(0,nrow = 32,ncol=8))
# set up data frame
dfRT <- data.frame(c(replicate(2, c("Fully Functional")), 
                     replicate(3, c("d = 1")), replicate(3, c("d = 2")), replicate(3, c("d = 3")),
                     replicate(3, c("d = 4")), replicate(3, c("d = 5")), replicate(3, c("d = 6")),
                     replicate(3, c("d = 7")), replicate(3, c("d = 8")), replicate(3, c("d = 9")),
                     replicate(3, c("d = 10"))),
                 c(c('kstar','Date'),replicate(10,c('ktilde','Date','TVE'))),
                 results)

# only needed if first column consists of numbers
dfRT[[1]] <- as.character(dfRT[[1]])

rle.lengths <- rle(dfRT[[1]])$lengths
first <- !duplicated(dfRT[[1]])
dfRT[[1]][!first] <- ""

# define appearance of \multirow
dfRT[[1]][first] <-
  paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", dfRT[[1]][first], "}}")

strCaption <- paste0("\\textbf{Simulation Results: } Table Caption goes here")

# set up xtable output
print(xtable(dfRT, digits = c(0, 0, 0, rep(0,8)), # first zero "represents" row numbers which we skip later
             align = "ll|c|cccccccc",  # align and put a vertical line (first "l" again represents column of row numbers)
             caption = strCaption, label = "SimTable"),
      size = "footnotesize", #Change size; useful for bigger tables "normalsize" "footnotesize"
      include.rownames = FALSE, #Don't print rownames
      include.colnames = FALSE, #We create them ourselves
      caption.placement = "top", #"top", NULL
      hline.after=NULL, #We don't need hline; we use booktabs
      floating=TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
      sanitize.text.function = force, # Important to treat content of first column as latex function
      add.to.row = list(pos = list(-1,
                                   nrow(dfRT)),
                        command = c(paste("\\toprule \n",  # NEW row
                                          #"\\multicolumn{3}{c}{}&\\multicolumn{5}{c}{\\textbf{Rand Index}} \\\\\n",
                                          #"\\cmidrule(l){4-8}\n",
                                          "Method & & Data 1 & Data 2 & Data 3 & Data 4 & Data 5 & Data 6 & Data 7 & Data 8 \\\\\n", # NEW row 
                                          "\\midrule \n"
                        ),
                        paste("\\bottomrule \n"  # paste is used as it is more flexible regarding adding lines
                        )
                        )
      )
)