# analyses done by long lin. just rerun his script with small modification because ngsremix is updated

#/home/genis/software/NGSremix/src/NGSremix -beagle /home/krishang/projects/miniprojects/maf_beagle/newbeagle.gz -fname /home/users/xi/Waterbuck_project/0.ColleaguesContribution/LinLong_Genis/RM_relate5Inds.waterbuck.DefassaWaterbuckRef.10.11.fopt_conv.gz -qname /home/users/xi/Waterbuck_project/0.ColleaguesContribution/LinLong_Genis/RM_relate5Inds.waterbuck.DefassaWaterbuckRef.10.11.qopt_conv -parental 1 -o tt -select 6,9,28,30,33,65,70,87,91,112,118

#/home/genis/software/NGSremix/src/NGSremix -beagle /home/krishang/projects/miniprojects/maf_beagle/newbeagle.gz -fname /home/users/xi/Waterbuck_project/0.ColleaguesContribution/LinLong_Genis/RM_relate5Inds.waterbuck.DefassaWaterbuckRef.10.11.fopt_conv.gz -qname /home/users/xi/Waterbuck_project/0.ColleaguesContribution/LinLong_Genis/RM_relate5Inds.waterbuck.DefassaWaterbuckRef.10.11.qopt_conv -bothanc 1 -o tt -select 6,9,28,30,33,65,70,87,91,112,118 -P 40
 
#/home/genis/github/NGSremix/src/NGSremix -beagle /home/krishang/projects/miniprojects/maf_beagle/newbeagle.gz -fname /home/users/xi/Waterbuck_project/0.ColleaguesContribution/LinLong_Genis/RM_relate5Inds.waterbuck.DefassaWaterbuckRef.10.11.fopt_conv.gz -qname /home/users/xi/Waterbuck_project/0.ColleaguesContribution/LinLong_Genis/RM_relate5Inds.waterbuck.DefassaWaterbuckRef.10.11.qopt_conv -bothanc 1 -o tt_2boot -select 6,9,28,30,33,65,70,87,91,112,118 -P 20 -boot 100 -block 1000


Rscript /home/genis/github/apoh/apoh.R -i /home/genis/africa1kg/waterbuck/recent_hybrids/run2/tt_2boot.anccoef -o run2_withboot --ids /home/genis/africa1kg/waterbuck/recent_hybrids/run2/waterbuck.ids


