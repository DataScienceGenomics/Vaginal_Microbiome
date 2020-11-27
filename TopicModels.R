# Extract observed counts
deng.counts <- as.data.frame(OTU_Table)
RemoveSparseFeatures(deng.counts)
##Is everything numeric
which(apply(deng.counts, 2, var) == 0)

#apply(deng.counts, 2, is.numeric)
Vag_Data_V34.FitGoM <- FitGoM(deng.counts, K=4, options="BF", tol=100, num_trials=5,
                              control=list(tmax=100))
# extract the omega matrix: membership weights of each cell
names(Vag_Data_V34.FitGoM$fit)
omega <- Vag_Data_V34.FitGoM$fit$omega
theta <- Vag_Data_V34.FitGoM$fit$theta

####Plot
#OmegaPlot
StructureGGplot(omega= omega,
                #annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Set1"),
                yaxis_label = "Vaginal Samples",
                order_sample = TRUE, axis_tick = list(axis_ticks_length = .5,
                                                      # axis_ticks_lwd_y = .1,
                                                      axis_ticks_lwd_x = .1,
                                                      axis_label_size = 7,
                                                      axis_label_face = "bold"))


##ThetaPlot
theta <- Vag_Data_V34.FitGoM$fit$theta
z<-as.data.frame(theta)

z$TaxaID<-rownames(z)
meltOTU_Table<-fast_melt1(OTU_T_glom)

ljthetaSpecies<-left_join(z,meltOTU_Table,by="TaxaID")
subsetSThetea<-ljthetaSpecies%>%dplyr::select(Species,`1`,`2`,`3`,`4`)
uniquesubsetSThetea<-unique(subsetSThetea)
rownames(uniquesubsetSThetea)<-uniquesubsetSThetea$Species
finalthetaSpecies<-uniquesubsetSThetea %>% dplyr::select(-Species)
pheatmap(finalthetaSpecies, scale="row", clustering_distance_rows="euclidean",
         fontsize_col=17,fontsize_row = 12,
)
graphics.off()
pdf("~/Desktop/SpeciesNamedHeatmap.pdf", width=10, height=18)
pheatmap(finalthetaSpecies, scale="row", clustering_distance_rows="euclidean",
         fontsize_col=17,fontsize_row = 12,
)
dev.off()
