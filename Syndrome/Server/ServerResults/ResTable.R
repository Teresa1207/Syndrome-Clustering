setwd("~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ServerResults")
source("Results_Serv.R")
require(xtable)

# set up data frame
dfRT <- data.frame(c(replicate(4, c("Category 1")), replicate(4, c("Category 2")),
                   replicate(4, c("Category 3")), replicate(2, c("Category 4")),
                   replicate(4, c("Category 5")), replicate(3, c("Category 6"))),
                 c(1:21),
                 ResTable[,-c(2:4)])

# only needed if first column consists of numbers
dfRT[[1]] <- as.character(dfRT[[1]])

rle.lengths <- rle(dfRT[[1]])$lengths
first <- !duplicated(dfRT[[1]])
dfRT[[1]][!first] <- ""

# define appearance of \multirow
dfRT[[1]][first] <-
  paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", dfRT[[1]][first], "}}")

strCaption <- paste0("\\textbf{Simulation Results: }For each scenario we calculated (from left) the percentage of times the correct number of clusters (3) was chosen using the CH index. The average Rand Index when K = 3 for 500 replications of our proposed syndrome based (SB)clustering method, compared to another one multidimensional method (KLM3d) and a univariate method using a single marker ($M_1,M_2,M_3$)")

# set up xtable output
print(xtable(dfRT, digits = c(0, 0, 0, rep(2,6)), # first zero "represents" row numbers which we skip later
             align = "ll|c|cccccc",  # align and put a vertical line (first "l" again represents column of row numbers)
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
                                          "\\multicolumn{3}{c}{}&\\multicolumn{5}{c}{\\textbf{Rand Index}} \\\\\n",
                                          "\\cmidrule(l){4-8}\n",
                                          " &$\\#$ & K=3 & SB & KLM3d & $M_1$ & $M_2$ & $M_3$ \\\\\n", # NEW row 
                                          "\\midrule \n"
                        ),
                        paste("\\bottomrule \n"  # paste is used as it is more flexible regarding adding lines
                        )
                        )
      )
)
