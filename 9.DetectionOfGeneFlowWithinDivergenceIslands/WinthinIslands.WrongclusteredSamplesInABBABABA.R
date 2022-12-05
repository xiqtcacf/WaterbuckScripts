######highlighting wrong clustered samples in Bar plot in ABBABABA plots
library(dplyr)
library(ggplot2)
setwd('/home/xi/Desktop/Fst_NJtree_PCAngsd_Heter/outliers_ABBABABA')
###ccd Nairobi434 as an example
Bohor_reedbuck <- read.table("flip.ccd.Fixed.name.Goat.108inds.8outgroups.Bohor_reedbuck.557.list.abbababa.txt",header=T)
tt <- Bohor_reedbuck %>% group_by(H1_Pops.H2_Pops.H3_Pops)
highlight_tt1 <- Bohor_reedbuck %>% 
  filter(Bohor_reedbuck$Z > 3  |  Bohor_reedbuck$Z < -3 )
tt1 <- highlight_tt1 %>% group_by(H1_Pops.H2_Pops.H3_Pops)

ccd_401_1 <- read.table("434.ccd",header=T)
ccd_401 <- subset(ccd_401_1, Z > 3 | Z < -3)
ccd_401 <- ccd_401 %>% group_by(H1_Pops.H2_Pops.H3_Pops)

png(filename = "Nairobi434H2.png", width = 700,height=700,)
ggplot(Bohor_reedbuck,aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat)) +
  geom_boxplot(width = 0.4, colour="black")  +
  geom_point(data=tt1, 
             aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= '|Zscore| > 3'), 
             color='dark red',
             size=0.01) + guides(size=guide_legend("Source", override.aes=list(shape=15, size = 10))) +
  geom_point(data=ccd_401, 
             aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= 'WrongClusterSample 401 in PCAngsd'), 
             color='light green',
             size=1) + guides(size=guide_legend("Source", override.aes=list(shape=15, size = 10))) +
  #  geom_text(data=ccd_401,
  #            aes(y=reorder(H1_Pops.H2_Pops.H3_Pops, Dstat),x=Dstat, 
  #                label=ccd_401$H1_PopsInd.H2_PopsInd.H3_PopsInd, hjust=1.65, vjust=1, angle=-30, check_overlap = T)) +
  labs(y='Common-Common-Defassa-BohorReedbuck',x='D statistic') + 
  geom_vline(xintercept = 0, linetype = "dashed") + xlim(-0.4,0.4) + theme_light() +
  theme(legend.position = "top")+  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ 
  #  geom_point(data = tt1, aes(size="|Zscore| > 3", shape = NA), colour = "dark blue") + 
  geom_point(data = ccd_401, aes(size="WrongClusterSample 434(NairobiH2) in PCAngsd", shape = NA), colour = "light green") +
  guides(size=guide_legend("", override.aes=list(shape=16, size = 3))) + theme(legend.position = c(0.23,0.90))
dev.off()

####ddc  Ugalla5618 as an example
Bohor_reedbuck <- read.table("flip.ddc.Fixed.name.Goat.108inds.8outgroups.Bohor_reedbuck.557.list.abbababa.txt",header=T)

tt <- Bohor_reedbuck %>% group_by(H1_Pops.H2_Pops.H3_Pops)

highlight_tt1 <- Bohor_reedbuck %>% 
  filter(Bohor_reedbuck$Z > 3  |  Bohor_reedbuckto ask $Z < -3 )
tt1 <- highlight_tt1 %>% group_by(H1_Pops.H2_Pops.H3_Pops)

ddc_401_1 <- read.table("5618.ddc",header=T)
ddc_401 <- subset(ddc_401_1, Z > 3 | Z < -3)
ddc_401 <- ddc_401 %>% group_by(H1_Pops.H2_Pops.H3_Pops)

png(filename = "Ugalla5618H2.png", width = 700,height=700,)
ggplot(Bohor_reedbuck,aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat)) +
  geom_boxplot(width = 0.4, colour="black")  +
  geom_point(data=tt1, 
             aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= '|Zscore| > 3'), 
             color='dark red',
             size=0.01) + guides(size=guide_legend("Source", override.aes=list(shape=15, size = 10))) +
  geom_point(data=ddc_401, 
             aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= 'outliers in PCAngsd'), 
             color='light green',
             size=1) + guides(size=guide_legend("Source", override.aes=list(shape=15, size = 10))) +
  #  geom_text(data=ddc_401,
  #            aes(y=reorder(H1_Pops.H2_Pops.H3_Pops, Dstat),x=Dstat, 
  #                label=ddc_401$H1_PopsInd.H2_PopsInd.H3_PopsInd, hjust=1.1, vjust=1, angle=-30, check_overlap = T)) +
  labs(y='Defassa-Defassa-Common-BohorReedbuck',x='D statistic') + 
  geom_vline(xintercept = 0, linetype = "dashed") + xlim(-0.4,0.4) + theme_light() +
  theme(legend.position = "top")+  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ 
  #  geom_point(data = tt1, aes(size="|Zscore| > 3", shape = NA), colour = "dark blue") + 
  geom_point(data = ddc_401, aes(size="WrongClusterSample 5618(UgallaH2) in PCAngsd", shape = NA), colour = "light green") +
  guides(size=guide_legend("", override.aes=list(shape=16, size = 3))) + theme(legend.position = c(0.23,0.90))
dev.off()