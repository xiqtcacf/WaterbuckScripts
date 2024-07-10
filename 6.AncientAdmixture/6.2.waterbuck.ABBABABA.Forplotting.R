#example for plotting ABBABABA including 11 recently admixed samples, ABBABABA excluding recently admixed samples same but different input files.

########get the right format of abbababa results
setwd('/home/xi/Desktop/SysBio-waterbuck-Revision1-August-2023/revised_analyses/ABBABABA_11admixS')
#df.1 <- read.table("8outgroup_crashed_duplicate_within_samples_Waterbuck.species.pops.txt",header=T)
#df.2 <- read.table("H1",header=T)
#df.2$id  <- 1:nrow(df.2)
#out  <- merge(df.2,df.1, by = "Inds")
#out1 <- out[order(out$id), ]
#write.table(out1, "H1.pops.species", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)

###ccd
library(ggplot2)
library(dplyr)
Bohor_reedbuck <- read.table("flip.ccd.Fixed.name.Goat.108inds.8outgroups.Bohor_reedbuck.557.list.abbababa.txt",header=T)
Bohor_reedbuck_adm <- read.table("FinalUsing.ccd.all",header=T)
hist(Bohor_reedbuck$Z)
Bohor_reedbuck$p <- 2*pnorm(-abs(Bohor_reedbuck$Z))
Bohor_reedbuck$bonferroni_p =
  p.adjust(Bohor_reedbuck$p,
           method = "bonferroni")
tt <- Bohor_reedbuck %>% group_by(H1_Pops.H2_Pops.H3_Pops)
highlight_tt1 <- Bohor_reedbuck %>% 
  filter(Bohor_reedbuck$bonferroni_p < 0.05 )
tt1 <- highlight_tt1 %>% group_by(H1_Pops.H2_Pops.H3_Pops)
tt2 <- Bohor_reedbuck_adm %>% group_by(H1_Pops.H2_Pops.H3_Pops)
p <- 
  ggplot(Bohor_reedbuck,aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat)) +
  geom_boxplot(width = 0.4, colour="black")  + 
  geom_point(data=tt1, 
             aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= "Bonferroni corrected P < 0.05"), 
             color="darkred",
             size=0.01) +
  geom_point(data=tt2, 
               aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= "Recently admixed 11 samples with other individuals"), 
               color='blue',
               size=0.01) +
  labs(y='Common-Common-Defassa-BohorReedbuck',x='D statistic') + 
  geom_vline(xintercept = 0, linetype = "dashed") + xlim(-0.4,0.4) + theme_light() +
  theme(legend.position = "top")+  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ 
  geom_point(data = tt1, aes(size="Bonferroni corrected P < 0.05", shape = NA), colour = "darkred")+
  geom_point(data = tt2, aes(size="Recently admixed 11 samples with other individuals", shape = NA), colour = "blue")+
  guides(size=guide_legend("", override.aes=list(shape=16, size = 3))) + theme(legend.position = c(0.23,0.94))
ggsave("revised.Highlighting11Admix.Goat.Common-Common-Defassa-BohorReedbuck.pdf", p, height = 10, width = 10)
ggsave("revised.Highlighting11Admix.Goat.Common-Common-Defassa-BohorReedbuck.png", p, height = 10, width = 10)

####ddc
Bohor_reedbuck <- read.table("flip.ddc.Fixed.name.Goat.108inds.8outgroups.Bohor_reedbuck.557.list.abbababa.txt",header=T)
Bohor_reedbuck_adm <- read.table("FinalUsing.ddc.all",header=T)
hist(Bohor_reedbuck$Z)
Bohor_reedbuck$p <- 2*pnorm(-abs(Bohor_reedbuck$Z))
Bohor_reedbuck$bonferroni_p =
  p.adjust(Bohor_reedbuck$p,
           method = "bonferroni")
tt <- Bohor_reedbuck %>% group_by(H1_Pops.H2_Pops.H3_Pops)
highlight_tt1 <- Bohor_reedbuck %>% 
  filter(Bohor_reedbuck$bonferroni_p < 0.05 )
tt1 <- highlight_tt1 %>% group_by(H1_Pops.H2_Pops.H3_Pops)
tt2 <- Bohor_reedbuck_adm %>% group_by(H1_Pops.H2_Pops.H3_Pops)
p1 <- 
  ggplot(Bohor_reedbuck,aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat)) +
  geom_boxplot(width = 0.4, colour="black")  + 
  geom_point(data=tt1, 
             aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= 'Bonferroni corrected P < 0.05'), 
             color='dark red',
             size=0.01) + guides(size=guide_legend("Source", override.aes=list(shape=15, size = 10))) +
    geom_point(data=tt2, 
               aes(y=reorder(H1_Pops.H2_Pops.H3_Pops,Dstat),x=Dstat, color= "Recently admixed 11 samples with other individuals"), 
               color='blue',
               size=0.01) +
  labs(y='Defassa-Defassa-Common-BohorReedbuck',x='D statistic') + 
  geom_vline(xintercept = 0, linetype = "dashed") + xlim(-0.4,0.4) + theme_light() +
  theme(legend.position = "left") + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ 
  geom_point(data = tt1, aes(size="Bonferroni corrected P < 0.05", shape = NA), colour = "dark red")+ 
  geom_point(data = tt2, aes(size="Recently admixed 11 samples with other individuals", shape = NA), colour = "blue")+ 
  guides(size=guide_legend("", override.aes=list(shape=16, size = 3))) + theme(legend.position = c(0.23,0.94))

ggsave("revised.Highlighting11Admix.Goat.Defassa-Defassa-Common-BohorReedbuck.pdf", p1, height = 10, width = 10)
ggsave("revised.Highlighting11Admix.Goat.Defassa-Defassa-Common-BohorReedbuck.png", p1, height = 10, width = 10)
