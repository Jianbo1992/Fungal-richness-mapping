# To build a structural equation model

#e.g., fungal saprotrophs

#library(pacman)
pacman::p_load(tidyverse, openxlsx)

data <- read.xlsx("F:/ SEM_data.xlsx")
head(data)

dat <- data %>% dplyr::select(Fungal.saprotrophs, latitude, vegetation.classification, gpp, 
                              band1,band3,band7,band15,band17, cf, thick,cec,ph,soc,tk,tp, slope, elevation,sw)
head(dat)
dat$vegetation.classification <- as.factor(dat$vegetation.classification)
str(dat)
nrow(dat)

######
# Step 1: Select the driving forces
pacman::p_load(randomForest, rfUtilities, rfPermute, MuMIn)

#####################################
########## 1.1 Selection of climatic variables

climdat <- dat %>% dplyr::select(Fungal.saprotrophs, band1, band3,band7,band15,band17)
set.seed(123)
clim_rf <- randomForest(Fungal.saprotrophs ~ .,data=climdat,importance=TRUE,proximity=TRUE)
clim_rf

#Climatic variables significance by random forest model
set.seed(123)
clim_perm <- rf.significance(clim_rf, climdat, nperm=99, ntree=500)
clim_perm

#To check the importance of each climatic variable and map the importance of climatic variables
set.seed(123)
clim_rfP <- rfPermute(Fungal.saprotrophs ~ ., data = climdat, ntree = 500,
                      na.action = na.omit, nrep = 100, num.cores = 1)
clim_imp <- rfPermute::importance(clim_rfP, sort.by = NULL, decreasing = TRUE)

important_clim <- clim_imp %>% as_tibble(rownames = "names") %>%
  data.frame() %>% mutate(
    label = if_else(X.IncMSE.pval < 0.001,"***",
                    if_else(X.IncMSE.pval <0.01,"**",
                            if_else(X.IncMSE.pval<0.05,"*","ns"))),
    X.IncMSE = as.numeric(X.IncMSE)) %>% arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))

# Mapping the importance of climatic variables
ggplot(important_clim, aes(x = names, y = X.IncMSE)) +
  geom_bar(aes(fill = group),stat = "identity") +
  geom_text(aes(y = X.IncMSE + 1, label = label)) +
  labs(x = "", y = "%IncMSE") + coord_flip() + theme_bw()

# To select the climatic variables that significant affect the richness of fungal saprotrophs
sig_clim <- important_clim %>% filter(group == "Sig") %>% .$names
finalclim <- climdat %>% dplyr::select(as.character(sig_clim))
finalclim


#####
########## 1.2 Selection of soil variables

soildat <- dat %>% dplyr::select(Fungal.saprotrophs, cf, thick,cec,ph,soc,tk,tp,sw)
set.seed(123)
soil_rf <- randomForest(Fungal.saprotrophs ~ .,data=soildat,importance=TRUE,proximity=TRUE)
soil_rf

# Check the significance of soil variables by random forest model
set.seed(123)
soil_perm <- rf.significance(soil_rf, soildat, nperm=99, ntree=500)
soil_perm

# To check the importance of each soil variable and map the importance of soil variables
set.seed(123)
soil_rfP <- rfPermute(Fungal.saprotrophs ~ ., data = soildat, ntree = 500,
                      na.action = na.omit, nrep = 100, num.cores = 1)
soil_imp <- rfPermute::importance(soil_rfP, sort.by = NULL, decreasing = TRUE)

important_soil <- soil_imp %>% as_tibble(rownames = "names") %>%
  data.frame() %>% mutate(
    label = if_else(X.IncMSE.pval < 0.001,"***",
                    if_else(X.IncMSE.pval <0.01,"**",
                            if_else(X.IncMSE.pval<0.05,"*","ns"))),
    X.IncMSE = as.numeric(X.IncMSE)) %>% arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))

# Mapping the importance of soil variables
ggplot(important_soil, aes(x = names, y = X.IncMSE)) +
  geom_bar(aes(fill = group),stat = "identity") +
  geom_text(aes(y = X.IncMSE + 1, label = label)) +
  labs(x = "", y = "%IncMSE") + coord_flip() + theme_bw()

# To select the soil variables that significant affect the richness of fungal saprotrophs
sig_soil <- important_soil %>% filter(group == "Sig") %>% .$names

#finalsoil <- soildat %>% dplyr::select(as.character(sig_soil))
finalsoil <- soildat %>% dplyr::select(ph,tp)


########## 1.3 Selection of topographic variables

topodat <- dat %>% dplyr::select(Fungal.saprotrophs, slope, elevation)
set.seed(123)
topo_rf <- randomForest(Fungal.saprotrophs ~ .,data=topodat,importance=TRUE,proximity=TRUE)
topo_rf

# Check the importance of each topographic variable
set.seed(123)
topo_perm <- rf.significance(topo_rf, topodat, nperm=99, ntree=500)
topo_perm
# Due to the results of random forest model are not significant (p = 0.96 and Var explained: -19.12), so we used correlation analysis to select topographic variables.

#library(BiocManager)
#BiocManager::install("dredge")
#library(dredge)
topolm <- lm(Fungal.saprotrophs ~ ., topodat)
options(na.action = "na.fail")
mod <- dredge
topolm
summary(topolm)

##########vegetation variable
finalveg <- vegdat %>% dplyr::select(gpp)

###### 
# Step 2. Integrated various variables based on principal component analysis
pacman::p_load(FactoMineR, factoextra)
# 2.1 climatic variables
finalclim
climpca <- PCA(finalclim, scale.unit = TRUE, graph = FALSE)

get_eigenvalue(climpca)
get_pca_var(climpca)$contrib
get_pca_var(climpca)$cor

fviz_pca_var(climpca, col.var = "contrib", repel = TRUE, 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# The PCA result of climatic variables
climind <- get_pca_ind(climpca)
climPC <- climind$coord[, 1:3]
colnames(climPC) <- c("climPC1", "climPC2", "climPC3")
head(climPC)
summary(climPC)

# 2.2 soil variables
finalsoil
soilpca <- PCA(finalsoil, scale.unit = TRUE, graph = FALSE)
get_eigenvalue(soilpca)
get_pca_var(soilpca)$contrib
get_pca_var(soilpca)$cor
fviz_pca_var(soilpca, col.var = "contrib", repel = TRUE, 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# The PCA result of soil variables
soilind <- get_pca_ind(soilpca)
soilPC <- soilind$coord[, 1, drop = F]
colnames(soilPC) <- "soilPC1"
head(soilPC)

# 2.2 topographic variables
finaltopo <- topodat[, -1]
topopca <- PCA(finaltopo, scale.unit = TRUE, graph = FALSE)
get_eigenvalue(topopca)
get_pca_var(topopca)$contrib
get_pca_var(topopca)$cor
fviz_pca_var(topopca, col.var = "contrib", repel = TRUE, 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# The PCA result of topographic variables
topoind <- get_pca_ind(topopca)
topoPC <- topoind$coord[, 1, drop = F]
colnames(topoPC) <- "topoPC1"
head(topoPC)


finaldat <- cbind.data.frame(finalclim, finalsoil, finaltopo, finalveg)
PCdat <- cbind.data.frame(climPC, soilPC, topoPC)

# Step 3. To build a SEM
pacman::p_load(piecewiseSEM, tidyverse, emmeans)
finaldat <- cbind.data.frame(richness = dat$Fungal.saprotrophs, latitude = dat$latitude, 
                             climPC, soilPC, topoPC, veg = finalveg)
names(finaldat)

fit <- psem(lm(richness ~ gpp + climPC2 + soilPC1, data = finaldat), 
            lm(soilPC1 ~ gpp + climPC2 + latitude, data = finaldat),
            lm(gpp ~ climPC2 + latitude, data = finaldat),
            lm(climPC2 ~ latitude+topoPC1, data = finaldat))
summary(fit)

