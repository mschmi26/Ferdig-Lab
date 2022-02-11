#PCA adapted from Sthda Principal Component Analysis in R: prcomp vs princomp
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
library(dplyr)
setwd('/Users/macke/Documents/PCA')
eastman  = read.csv("drugBatchMap_PCA.csv", 
                    header = TRUE, 
                    row.names = 1, 
                    stringsAsFactors = FALSE, 
                    fileEncoding="UTF-8-BOM")

eastman_NOna= eastman[ , colSums(is.na(eastman)) == 0]
#install.packages("factoextra")
library(factoextra)
### load in data and 
#library("factoextra")
data(decathlon2)
decathlon2.active <- decathlon2[1:23, 1:10]
head(decathlon2.active[, 1:6])

decathlon2.active <- decathlon2[1:23, 1:10]
decathlon2.active = eastman_NOna

sapply(eastman, function(x) sum(is.na(x)))
rowSums(is.na(eastman))


### compute PCA and scree plot
res.pca <- prcomp(decathlon2.active, scale = FALSE)
fviz_screeplot(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

### graph variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
### biplot individuals and variables
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

###qualtitative / categortical variables
### grouping individuals
groups <- as.factor(decathlon2$Competition[1:23])
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)


library(magrittr) # for pipe %>%
library(dplyr)   # everything else
# 1. Individual coordinates
res.ind <- get_pca_ind(res.pca)
# 2. Coordinate of groups
coord.groups <- res.ind$coord %>%
  as_tibble() %>%
  select(Dim.1, Dim.2) %>%
  mutate(competition = groups) %>%
  group_by(competition) %>%
  summarise(
    Dim.1 = mean(Dim.1),
    Dim.2 = mean(Dim.2)
  )
coord.groups


###Theory begind PCA



